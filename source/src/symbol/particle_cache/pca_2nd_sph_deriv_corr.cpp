/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018,
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
 * and others authors stated in the AUTHORS file in the top-level 
 * source directory.
 *
 * SYMPLER is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SYMPLER is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SYMPLER.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Please cite the research papers on SYMPLER in your own publications. 
 * Check out the PUBLICATIONS file in the top-level source directory.
 *
 * You are very welcome to contribute extensions to the code. Please do 
 * so by making a pull request on https://github.com/kauzlari/sympler
 * 
 */


#include "pca_2nd_sph_deriv_corr.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "math_helper.h"
#include "gsl_helper.h"
#include "threads.h"


const SymbolRegister<PCa2ndSPHDerivCorr> pca_2nd_sph_deriv_corr("SPH2ndDerivCorr");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

// Macro used to avoid code repetition. myOrOther = my increments the my*
// temporary arrays, myOrOther = other the other* arrays. Note that the
// "string" replacements are done at preprocessor time, i.e., far before
// runtime. Meaningful argument combinations so far: (my, 1.0) or (other, -1.0)
#define INCREMENT_ARRAYS(myOrOther, macroAntiSymmetryFactor) \
{\
						totalAntiSymmetryFactor \
							= outerAntiSymmetryFactor * macroAntiSymmetryFactor; \
						/* contribution to matrix E step 1 ( / distSq missing!) */ \
						tempEabe \
							= totalAntiSymmetryFactor * \
							myOrOther##MinVolumejTimesWeight * baseProd; \
\
						/* contribution to matrix G^T*/ \
						/* see docu above for the following identity */ \
						(*myOrOther##T_Gcab_abc)[helperIdxG] += tempEabe; \
\
						/* contribution to matrix E step 2: / distSq */ \
						tempEabe /= distSq; \
						(*myOrOther##Eabe)[helperIdxG] += tempEabe; \
\
						/* Caacc terms: row = (a, b), column = c */ \
						/* Yes we can also use helperIdx here! */ \
						/* Caacc = Eaac x distVecc */ \
						(*myOrOther##Cabcd)[helperIdxCcc] \
						/* Negate effect of totalAntiSymmetryFactor for the symmetric \
						 * Cabcd */ \
							+= totalAntiSymmetryFactor \
								* tempEabe * distVec[c]; \
\
						/* d loop for the half of the c != d terms which we */ \
						/* compute */ \
						for (size_t d = c+1; d < SPACE_DIMS; ++d) { \
\
							/* a determines the row (1 of the first 3) */ \
							/*NO!!!! It could also be an ab-row !!!;*/ \
							/* for the column we have helperColumnIdxC */ \
							/*NO, we try to take helperIdxCcc which seems to be valid no */ \
							/*matter if in an aa or ab case. We just subtract c and add */ \
							/*myOrOther##Helper_dDeltaIdx for c != d */ \
							(*myOrOther##Cabcd) \
							/*[myOrOther##HelperColumnIdxC + aRowCompleteIdxC]*/ \
							[myOrOther##Helper_dDeltaIdx - c + helperIdxCcc] \
								+= totalAntiSymmetryFactor \
									* tempEabe * distVec[d]; \
\
							/*++myOrOther##HelperColumnIdxC;*/ \
							++myOrOther##Helper_dDeltaIdx; \
\
						} /* end of for (size_t d... */ \
} \
while(0)


const size_t PCa2ndSPHDerivCorr::s_matEntries
	= SPACE_DIMS * (SPACE_DIMS + 1) / 2;

// definition of and memory allocation for static member variable
gsl_vector* PCa2ndSPHDerivCorr::s_rhs
	= gsl_vector_alloc(PCa2ndSPHDerivCorr::s_matEntries);

// Define our static initialiser, which will call the _Initialiser() default
// constructor, which will initialize PCa2ndSPHDerivCorr::s_rhs
PCa2ndSPHDerivCorr::_Initialiser PCa2ndSPHDerivCorr::s_initialiser;



PCa2ndSPHDerivCorr::PCa2ndSPHDerivCorr(/*Node*/Simulation* parent)
  : ParticleCache(parent)
{
  m_datatype = DataFormat::TENSOR;

	m_workingMat.size1 = m_workingMat.size2 = m_workingMat.tda = s_matEntries;
  m_workingMat.owner = 0;
  // will be assigned in computeCacheFor
  m_workingMat.data = NULL;

  m_permutation = gsl_permutation_alloc(s_matEntries);

  m_B = gsl_vector_alloc(s_matEntries);

  init(); 
}

PCa2ndSPHDerivCorr::~PCa2ndSPHDerivCorr()
{
  gsl_vector_free(m_B);
  gsl_permutation_free(m_permutation);
}

void PCa2ndSPHDerivCorr::init()
{

  m_properties.setClassName("PCa2ndSPHDerivCorr");

  m_properties.setName("PCaSPH2ndDerivCorr");

  m_properties.setDescription
	(
		"Computes the symmetric 3x3 correction matrix B for the second SPH-"
		"derivative according to R. Fatehi, M.T. Manzari, Computers and "
		"Mathematics with Applications 61 (2011) 482–498. B is obtained from the "
		"equation B:A == -I, where A is a 3x3x3x3 matrix. I is the 3x3 identity matrix."
  );
      
  STRINGPC (weightingFunction, m_wfName,
  	"Symbol name of the externally defined weighting function to be used. "
    "SPH2ndDerivCorr will use the derivative of the defined weighting "
  	"(interpolation) function. Also note that the pair-computations within "
  	"this module will use the cutoff distance defined by the used weighting "
  	"function.");

  m_wfName = "undefined";

  STRINGPC (volume, m_volumeName,
  	"Symbol name of the externally computed SPH particle volume.");

  m_volumeName = "undefined";

  STRINGPC (SPH1stDerivCorr, m_1stDerivCorrName,
  	"Symbol name of the externally computed correction matrix for the SPH-"
    "discretisation of the 1st derivative.");

  m_1stDerivCorrName = "undefined";

  STRINGPC (symbol, m_symbolName,
  	"Symbol name of the correction matrix for the SPH-discretisation of the "
  	"2nd derivative computed by this module.");

  m_symbolName = "B2";

  // Overriding description to make clear that we do not allow for
  m_properties.setPropDescription
	("species",
  		"Name for the species of the particles, this Symbol is used for. "
  		"'species = \"ALL\"' is not allowed for this module."

  				);
}

void PCa2ndSPHDerivCorr::setup()
{

	// we skip the direct parent's ParticleCache::setup() since some of the
	// parent's features are switched off here, such as m_species = "ALL"
	Symbol::setup();

  // turns off the default gsl_error handling so that we can handle ourselves
  gsl_set_error_handler_off ();

  if(m_species == "undefined")
    throw gError("PCa2ndSPHDerivCorr::setup", "Attribute 'species' has value "
    		"\"undefined\".");

  if(m_symbolName == "undefined")
    throw gError("PCa2ndSPHDerivCorr::setup", "Attribute 'symbol' has value "
    		"\"undefined\".");

  if(m_wfName == "undefined")
    throw gError("PCa2ndSPHDerivCorr::setup", "Attribute 'weightingFunction' "
    		"has value \"undefined\".");

  if(m_volumeName == "undefined")
    throw gError("PCa2ndSPHDerivCorr::setup", "Attribute 'volume' "
    		"has value \"undefined\".");

  if(m_1stDerivCorrName == "undefined")
    throw gError("PCa2ndSPHDerivCorr::setup", "Attribute 'SPH1stDerivCorr' "
    		"has value \"undefined\".");

  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw
			gError
				("PCa2ndSPHDerivCorr::setup", "Attribute 'stage' has none of the "
						"allowed values \"0\", \"1\", \"2\".");

  
  // should we create a Cache for the other colours too?
  if(m_species == "ALL") {

  	throw gError(
  			"PCa2ndSPHDerivCorr::setup", "'species = \"ALL\"' is not allowed in "
				"this module. Please give one specific species and repeat for other "
				"species if needed.");

  } // end of if(m_species == "ALL")

  else {
    m_colour = M_MANAGER->getColour(m_species);

    m_wf = M_SIMULATION->findWeightingFunction(m_wfName);

    checkInputSymbolExistences(m_colour);

    // is the symbol already existing?
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      throw gError("PCa2ndSPHDerivCorr::setup", "Symbol " + m_symbolName +
      		" is already existing for species '" + M_MANAGER->species(m_colour) +
					"'. Second definition is not allowed in PCa2ndSPHDerivCorr.");

    m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName,
    		m_datatype, false, m_symbolName).offset;

    // system matrix
    if(Particle::s_tag_format[m_colour].attrExists("__sysMatFor" + m_symbolName))
      throw gError("PCa2ndSPHDerivCorr::setup", "Symbol \"__sysMat\" is "
      	"already existing for species '" + M_MANAGER->species(m_colour) +
				"'. But this name is reserved for module SPH2ndDerivCorr and can only "
				"be used once. Aborting.");

    // !!! NOTE: All DataFormat::VECTOR_DOUBLE attributes MUST be persistent
    // (therefore the "true") because they will otherwise be cleared once per
    // time step, which is fatal in this case because it clears all the
    // contents, i.e., the whole vector<double> to which the vector_double_sp
    // is pointing to by even freeing the allocated memory and does not just
    // zeroing the entries. FIXME: Check if this is the most meaningful
    // behaviour or if we should weaken the clear() procedure for SmartPointers

  	m_systemMatOffset
			= Particle::s_tag_format[m_colour].addAttribute("__sysMatFor" + m_symbolName,
					DataFormat::VECTOR_DOUBLE, true, "__sysMat").offset;


    // the 1st derivative correction MUST already exist for m_colour
  	if(!Particle::s_tag_format[m_colour].attrExists(m_1stDerivCorrName))
  		throw gError("PCa2ndSPHDerivCorr::setup", "Symbol " + m_1stDerivCorrName +
  				" not found for species '" + M_MANAGER->species(m_colour) + "' but "
					"required within this module.");

  	// Does it have the required DataFormat?
  	if (Particle::s_tag_format[m_colour].attrByName(m_1stDerivCorrName).datatype
  			!= DataFormat::TENSOR)
  		throw gError("PCa2ndSPHDerivCorr::setup", "Symbol " + m_1stDerivCorrName +
  				" found for species '" + M_MANAGER->species(m_colour) + "' but "
					"with wrong data format \""
					+ Particle::s_tag_format[m_colour].attrByName(m_1stDerivCorrName)
					.datatypeAsString() + "\" different from a 3x3 tensor.");

    // set memory offset
  	m_1stDerivCorrOffset
			= Particle::s_tag_format[m_colour].attrByName(m_1stDerivCorrName)
			.offset;


  	// START: creating the helper matrices

    if(Particle::s_tag_format[m_colour].attrExists("__CabcdFor" + m_symbolName))
      throw gError("PCa2ndSPHDerivCorr::setup", "Symbol \"__Cabcd\" is "
      	"already existing for species '" + M_MANAGER->species(m_colour) +
				"'. But this name is reserved for module SPH2ndDerivCorr and can only "
				"be used once. Aborting.");

    m_CabcdHelperOffset
			= Particle::s_tag_format[m_colour].addAttribute("__CabcdFor" + m_symbolName,
					DataFormat::VECTOR_DOUBLE, true, "__Cabcd").offset;


    if(Particle::s_tag_format[m_colour].attrExists("__HabcFor" + m_symbolName))
      throw gError("PCa2ndSPHDerivCorr::setup", "Symbol \"__Habc\" is "
      	"already existing for species '" + M_MANAGER->species(m_colour) +
				"'. But this name is reserved for module SPH2ndDerivCorr and can only "
				"be used once. Aborting.");

    m_HabfHelperOffset
			= Particle::s_tag_format[m_colour].addAttribute("__HabcFor" + m_symbolName,
					DataFormat::VECTOR_DOUBLE, true, "__Habc").offset;


    if(Particle::s_tag_format[m_colour].attrExists("__EabcFor" + m_symbolName))
      throw gError("PCa2ndSPHDerivCorr::setup", "Symbol \"__Eabc\" is "
      	"already existing for species '" + M_MANAGER->species(m_colour) +
				"'. But this name is reserved for module SPH2ndDerivCorr and can only "
				"be used once. Aborting.");

    m_EabeHelperOffset
			= Particle::s_tag_format[m_colour].addAttribute("__EabcFor" + m_symbolName,
					DataFormat::VECTOR_DOUBLE, true, "__Eabc").offset;


    if(Particle::s_tag_format[m_colour].attrExists("__TGabcFor" + m_symbolName))
      throw gError("PCa2ndSPHDerivCorr::setup", "Symbol \"__TGabc\" is "
      	"already existing for species '" + M_MANAGER->species(m_colour) +
				"'. But this name is reserved for module SPH2ndDerivCorr and can only "
				"be used once. Aborting.");

    m_T_Gcab_abcHelperOffset
			= Particle::s_tag_format[m_colour].addAttribute("__TGabcFor" + m_symbolName,
					DataFormat::VECTOR_DOUBLE, true, "__TGabc").offset;

  	// END: creating the helper matrices


  	// START: resize all vector-doubles

    Particle::s_tag_format[m_colour].vectorDoubleByOffset(m_systemMatOffset)
    		-> resize(s_matEntries * s_matEntries);

    Particle::s_tag_format[m_colour].vectorDoubleByOffset(m_CabcdHelperOffset)
    		-> resize(s_matEntries * s_matEntries);

    Particle::s_tag_format[m_colour].vectorDoubleByOffset(m_HabfHelperOffset)
    		-> resize(s_matEntries * SPACE_DIMS);

    Particle::s_tag_format[m_colour].vectorDoubleByOffset(m_EabeHelperOffset)
    		-> resize(s_matEntries * SPACE_DIMS);

    Particle::s_tag_format[m_colour].vectorDoubleByOffset(m_T_Gcab_abcHelperOffset)
    		-> resize(s_matEntries * SPACE_DIMS);

   	// END: resize all vector_doubles

  } // end of else of if(m_species == "ALL")

  // register for precomputation (= external call of
  // PCa2ndSPHDerivCorr::precompute())
  if(m_phaseUser == 0) {
  	M_CONTROLLER->registerForPrecomputation_0(this);
    Particle::registerCache_0(this);
  }
  else if(m_phaseUser == 1) {
  	M_CONTROLLER->registerForPrecomputation(this);
    Particle::registerCache(this);
  }
  else if(m_phaseUser == 2) {

  	// register a copy
    ParticleCache* pc = copyMySelf();

    Particle::registerCache(pc);
    Particle::registerCache_0(this);
    M_CONTROLLER->registerForPrecomputation(pc);
    M_CONTROLLER->registerForPrecomputation_0(this);
  }
  else
  	throw
			gError("PCa2ndSPHDerivCorr::setup","attribute 'stage' has none of "
					"the values \"0\", \"1\", or \"2\".");
}

void PCa2ndSPHDerivCorr::checkInputSymbolExistences(size_t colour)
{
	m_volumeOffset.resize(M_PHASE -> nColours());

	// the volume MUST already exist for ALL colours
	for (size_t col = 0; col < M_PHASE -> nColours(); ++col) {
		if(!Particle::s_tag_format[col].attrExists(m_volumeName))
			throw gError("PCa2ndSPHDerivCorr::setup", "Symbol " + m_volumeName +
					" not found for species '" + M_MANAGER->species(col) + "' but "
					"required within this module.");

		// Does it have the required DataFormat?
		if (Particle::s_tag_format[col].attrByName(m_volumeName).datatype
				!= DataFormat::DOUBLE)
			throw gError("PCa2ndSPHDerivCorr::setup", "Symbol " + m_volumeName +
					" found for species '" + M_MANAGER->species(col) + "' but "
					"with wrong data format \""
					+ Particle::s_tag_format[col].attrByName(m_volumeName)
					.datatypeAsString() + "\" different from a scalar.");

		// set memory offset
		m_volumeOffset[col]
									 = Particle::s_tag_format[col].attrByName(m_volumeName)
									 .offset;

	} // end of for (size_t col ...
}


void PCa2ndSPHDerivCorr::precompute()
{
	m_pairLoopToDo = true;
}

void PCa2ndSPHDerivCorr::initSystemMatrices() {

	FOR_EACH_PARTICLE_C__PARALLEL(M_PHASE, m_colour, NULL,

		// (i->tag.vectorDoubleByOffset(m_CabcdHelperOffset)) is of type
		// vector_double_sp
		// *((i->tag.vectorDoubleByOffset(m_CabcdHelperOffset))) is of type
		// vector<double>
		vector<double>& Cabcd =
			*((i->tag.vectorDoubleByOffset(m_CabcdHelperOffset)));

		vector<double>& Eabe =
			*((i->tag.vectorDoubleByOffset(m_EabeHelperOffset)));

		vector<double>& T_Gcab_abc =
			*((i->tag.vectorDoubleByOffset(m_T_Gcab_abcHelperOffset)));

		for (size_t j = 0; j < SPACE_DIMS_SQUARED * (SPACE_DIMS + 1) * (SPACE_DIMS + 1) / 4; ++j)
			Cabcd[j] = 0.;

		for (size_t j = 0; j < SPACE_DIMS_SQUARED * (SPACE_DIMS + 1) / 2; ++j)
			Eabe[j] = T_Gcab_abc[j] = 0.;

		// Habf can be initialised on the fly

	);

}

void PCa2ndSPHDerivCorr::buildSystemMatrices() {

	point_t dummy = {{{ 0, 0, 0 }}};

	size_t helperIdx;

	Particle *myParticle, *otherParticle;

  double cutoff = m_wf->cutoff();

	size_t myColour, otherColour;

	FOR_EACH_COLOUR_PAIR(M_MANAGER,

		bool myColourIsFirst = cp->firstColour() == m_colour;
		bool myColourIsSecond = cp->secondColour() == m_colour;
		bool bothColours = myColourIsFirst && myColourIsSecond;
		double outerAntiSymmetryFactor = 1.0;
		double totalAntiSymmetryFactor;

		// we need the other particle even if it is not of m_colour (to access
		// its volume) !
		if (myColourIsFirst || myColourIsSecond) {

			FOR_EACH_PAIR__PARALLEL(PCa2ndSPHDerivCorr, cp,

				const double& dist = pair -> abs();

				if (dist < cutoff) {

				if (myColourIsFirst) {

					myParticle = pair -> firstPart();

					if (myColourIsSecond)

						otherParticle = pair -> secondPart();

				}
				else {

					outerAntiSymmetryFactor = -1.0;

					myParticle = pair -> secondPart();
					otherParticle = pair -> firstPart();

				}

				// in any case ...
				myColour = myParticle -> c;
				otherColour = otherParticle -> c;

				// matrix A = C + D with
				// Cabcd = -volumej * weight
				// 	* (distVeca x distVecb x distVecc x distVecd) / dist^2
				//	= Eabc x distVecD	= Gabe / dist^2 x distVecD
				// Dabcd = Eabe . Fef . Gfcd
				// Eabe = Gabe / dist^2 (Symmetry in all index pairs!)
				// Fef = firstDerivCorrMatrixef
				// Gfcd =
				//	-volumej * weight * (distVecf x distVecc x distVecd)
				// where x denotes an outer product

				const point_t& distVec = pair -> cartesian();

				const double& distSq = pair -> absSquare();

				double minWeight = m_wf -> weight(pair, dummy) * (-1.);

				double myMinVolumejTimesWeight = otherParticle
					-> tag.doubleByOffset(m_volumeOffset[otherColour]) * minWeight;

				// only thing we do if other Particle not modified
				double otherMinVolumejTimesWeight = 0.;

				// we need to use pointers here because we can not declare references
				// without initialisation. References will be introduced when needed
				// (i->tag.vectorDoubleByOffset(m_CabcdHelperOffset)) is of type
				// vector_double_sp
				// (i->tag.vectorDoubleByOffset(m_CabcdHelperOffset)).value() is of
				// type vector<double>*
				vector<double>* myCabcd =
					(myParticle->tag.vectorDoubleByOffset(m_CabcdHelperOffset))
							.value();
				vector<double>* myEabe =
					(myParticle->tag.vectorDoubleByOffset(m_EabeHelperOffset))
							.value();
				vector<double>* myT_Gcab_abc =
					(myParticle->tag.vectorDoubleByOffset(m_T_Gcab_abcHelperOffset))
							.value();

				vector<double>* otherCabcd;
				vector<double>* otherEabe;
				vector<double>* otherT_Gcab_abc;

				if (bothColours) {

					otherCabcd =
						(otherParticle->tag.vectorDoubleByOffset(m_CabcdHelperOffset))
								.value();
					otherEabe =
						(otherParticle->tag.vectorDoubleByOffset(m_EabeHelperOffset))
								.value();
					otherT_Gcab_abc =
						(otherParticle->tag.
								vectorDoubleByOffset(m_T_Gcab_abcHelperOffset)).value();

					otherMinVolumejTimesWeight = myParticle
						-> tag.doubleByOffset(m_volumeOffset[myColour]) * minWeight;

				}

				size_t helperRowIdx = SPACE_DIMS;

				////// START: build the matrices C, E, G

				// Eabe is treated symmetric in a, b, but full in e, so we
				// make a s_matEntries x SPACE_DIMS matrix out of it where the row represents the
				// two indices a, b in the typical order (3D)
				// 00, 11, 22, 01, 02, 12

				// For C, each row and column go through two indices in this
				// order (3D):
				// 00, 11, 22, 01, 02, 12
				// a loop for the row
				for (size_t a = 0; a < SPACE_DIMS; ++a) {

					// increment to column index for (c,d) with c != d dealt with in a loop
					// over d within INCREMENT_ARRAYS. Therefore it starts at SPACE_DIMS
					size_t myHelper_dDeltaIdx = SPACE_DIMS;
					size_t otherHelper_dDeltaIdx = SPACE_DIMS;

					// when row a is complete, the next index in matrix C is...
					size_t aRowCompleteIdxC = a * s_matEntries;

					// c loop for the column
					for (size_t c = 0; c < SPACE_DIMS; ++c) {

						// aac(c) terms: row = a, column = c
						size_t helperIdxCcc = c + aRowCompleteIdxC;
						size_t helperIdxG = c + a * SPACE_DIMS;

						double baseProd = distVec[a] * distVec[a] * distVec[c];
						// next is used in the following macro
						double tempEabe;

						// the calls of the array increments
						INCREMENT_ARRAYS(my, 1.0);
						if (bothColours) {
							INCREMENT_ARRAYS(other, -1.0);
						}
					} // end of for (size_t c...

					// b loop for the half of the a != b terms which we compute
					for (size_t b = a+1; b < SPACE_DIMS; ++b) {

						// increment to column index for (c,d) with c != d dealt with in a loop
						// over d within INCREMENT_ARRAYS. Therefore it starts at SPACE_DIMS
						size_t myHelper_dDeltaIdx = SPACE_DIMS;
						size_t otherHelper_dDeltaIdx = SPACE_DIMS;

						// when row (a,b) in matrix C is complete, the next index
						// is...
						size_t abRowCompleteIdxC = helperRowIdx * s_matEntries;

						// c loop for the column
						for (size_t c = 0; c < SPACE_DIMS; ++c) {

							// Eabc terms: row = (a,b), column = c(c)
							size_t helperIdxG = c + helperRowIdx * SPACE_DIMS;
							size_t helperIdxCcc = c + abRowCompleteIdxC;

							// next two are used in the following macro
							double baseProd = distVec[a] * distVec[b] * distVec[c];
							double tempEabe;

							// the calls of the array increments
							// only baseProd is different from first call
							INCREMENT_ARRAYS(my, 1.0);
							if (bothColours) {
								INCREMENT_ARRAYS(other, -1.0);
							}
						} // end of for (size_t c...

						++helperRowIdx;

					} // end of for (size_t b...

				} // end of for (size_t a...

				////// END: build the matrices C, E, G

				} // end of if (dist < cutoff)

			); // end of FOR_EACH_PAIR_PARALLEL

		} // end of if cp contains m_colour

	); // end of FOR_EACH_COLOUR_PAIR

#undef INCREMENT_ARRAYS

	////// START: compute Dabcd = Eabe . Fef . Gfcd

	FOR_EACH_PARTICLE_C__PARALLEL(M_PHASE, m_colour, NULL,

		// compute 6x3 matrix Habf = Eabe . Fef
		double* Fef = (i->tag.tensorByOffset(m_1stDerivCorrOffset)).tensor;

		// access the Eabe-matrix in GSL form. The RHS must become a double*
		// (p->tag.vectorDoubleByOffset(m_EabeHelperOffset)) is of type
		// vector_double_sp
		// *((p->tag.vectorDoubleByOffset(m_EabeHelperOffset))) is of type
		// vector<double>
		// (*((p->tag.vectorDoubleByOffset(m_EabeHelperOffset))))[0] is of
		// type double
		// &((*((p->tag.vectorDoubleByOffset(m_EabeHelperOffset))))[0]) is
		// of type double*
		double* Eabe
			= &((*((i->tag.vectorDoubleByOffset(m_EabeHelperOffset))))[0]);

		// dito for Habf
		double* Habf
			= &((*((i->tag.vectorDoubleByOffset(m_HabfHelperOffset))))[0]);

		MathHelper::matrixMatrixProd
			(s_matEntries, SPACE_DIMS, SPACE_DIMS, Eabe, Fef, Habf);

		// compute Dabcd = Eabe . Fef . Gfcd = Habf . Gfcd = Habf . Gfcd^T^T

		// access the Aabcd-matrix in GSL form. The RHS must become a double*
		// (p->tag.vectorDoubleByOffset(m_systemMatOffset)) is of type vector_double_sp
		// *((p->tag.vectorDoubleByOffset(m_systemMatOffset))) is of type vector<double>
		// (*((p->tag.vectorDoubleByOffset(m_systemMatOffset))))[0] is of type double
		// &((*((p->tag.vectorDoubleByOffset(m_systemMatOffset))))[0]) is of type double*
		double* Aabcd
			= &((*((i->tag.vectorDoubleByOffset(m_systemMatOffset))))[0]);

		// dito for T_Gcab_abc
		double* T_Gcab_abc
			= &((*((i->tag.vectorDoubleByOffset(m_T_Gcab_abcHelperOffset))))[0]);

		MathHelper::matrixMatrixTProd
			(s_matEntries, s_matEntries, SPACE_DIMS, Habf, T_Gcab_abc, Aabcd);


		////// END: compute Dabcd = Eabe . Fef . Gfcd


		// Now A <- C + A = C + D = C + E.F.G
		double* Cabcd
			= &((*((i->tag.vectorDoubleByOffset(m_CabcdHelperOffset))))[0]);

		for (size_t _i = 0; _i < s_matEntries * s_matEntries; ++_i)
			Aabcd[_i] += Cabcd[_i];

	); // end of FOR_EACH_PARTICLE_C__PARALLEL

} // end of void PCa2ndSPHDerivCorr::buildSystemMatrices()

void PCa2ndSPHDerivCorr::computeCacheFor(Particle* p) {

	// precompute the 6x6 B-matrix for each particle.
	// (FIXME: what a waste! Improve Neighbour list to avoid that!)

	// B-matrix computation is done once in a pair loop for ALL particles inside
	// the first call of computeCacheFor for the first particle
	if(m_pairLoopToDo) {

		initSystemMatrices();
		buildSystemMatrices();

		m_pairLoopToDo = false;

	} // end of if(m_pairLoopToDo)

	// Now solve the linear system

	// access the A-matrix in GSL form. The RHS must become a double*
	// (p->tag.vectorDoubleByOffset(m_systemMatOffset)) is of type vector_double_sp
	// *((p->tag.vectorDoubleByOffset(m_systemMatOffset))) is of type vector<double>
	// (*((p->tag.vectorDoubleByOffset(m_systemMatOffset))))[0] is of type double
	// &((*((p->tag.vectorDoubleByOffset(m_systemMatOffset))))[0]) is of type double*
	m_workingMat.data
		= &((*((p->tag.vectorDoubleByOffset(m_systemMatOffset))))[0]);

  // Compute the LU decomposition
  int signum;
  int error = gsl_linalg_LU_decomp(&m_workingMat, m_permutation, &signum);

  if(error)
  	throw
			gError(
					"PCa2ndSPHDerivCorr::computeCacheFor",
					"gsl_linalg_LU_decomp failed with return value "
					+ ObjToString(error)
					);

  // Solve the system for B
  error = gsl_linalg_LU_solve(&m_workingMat, m_permutation, s_rhs, m_B);

  if(error) {

  	MSG_DEBUG("PCa2ndSPHDerivCorr::computeCacheFor",
  			"ERROR for system matrix (p=" << p -> mySlot << "): ");

    for (size_t i = 0; i < s_matEntries; i++) {
    	cout << endl;
      for (size_t j = 0; j < s_matEntries; j++)
          cout << gsl_matrix_get (&m_workingMat, i, j) << " ";
    }

    cout << endl << " (might have been modified by GSL LU routine), "
				"particle: " << p->mySlot << ", colour: " << p->c << endl;

  	if(gsl_linalg_LU_det(&m_workingMat, signum) == 0.0) {

  		throw
				gError("PCa2ndSPHDerivCorr::computeCacheFor",
						"failed because matrix is singular.");
  	}
  	else {

  		throw
				gError("PCa2ndSPHDerivCorr::computeCacheFor",
						"gsl_linalg_LU_solve failed with error "
						+ string(gsl_strerror(error)));
  	}
  } // end of if(error)

  // Copy solution to tag-space

	tensor_t& Bmat = p->tag.tensorByOffset(m_offset);

	size_t idx = SPACE_DIMS;

  for (size_t i = 0; i < SPACE_DIMS; ++i) {
  	// diagonal elements
  	Bmat(i, i) = GSL_VECTOR_GET(m_B, i);

  	// off-diagonal elements
  	for (size_t j = i+1; j < SPACE_DIMS; ++j) {
      // 0.5 due to storage convention for gsl_vector B (see doc)
  		Bmat(i, j) = Bmat(j, i) = 0.5*GSL_VECTOR_GET(m_B, idx);
  		++idx;
  	}
  }

} // end of PCa2ndSPHDerivCorr::computeCacheFor(..)
