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


#include "pca_2nd_sph_deriv_corrTest.h"
#include "fake_pca_2nd_sph_deriv_corr.h"

#include "cppunit_helper.h"
#include "gsl_helper.h"

CPPUNIT_TEST_SUITE_REGISTRATION (PCa2ndSPHDerivCorrTest);


void PCa2ndSPHDerivCorrTest :: setUp (void)
{
	// will call PCa2ndSPHDerivCorrTest :: setUpParticleCache()
	ParticleCacheTest::setUp();
}

void PCa2ndSPHDerivCorrTest :: setUpParticleCache (void)
{
	m_particleCache = new FakePCa2ndSPHDerivCorr(m_simulation);
}

//void PCa2ndSPHDerivCorrTest :: tearDown (void)
//{
//	ParticleCache::tearDown();
//}

void PCa2ndSPHDerivCorrTest :: initTest (void)
{

	ParticleCacheTest::initTest();

	CPPUNIT_ASSERT_EQUAL
		(string("PCaSPH2ndDerivCorr"),
				m_particleCache -> returnProperties().name());

	CppunitHelper::testPropAttr(
		string("weightingFunction"),
		string(
		  "Symbol name of the externally defined weighting function to be used. "
		  "SPH2ndDerivCorr will use the derivative of the defined weighting "
		  "(interpolation) function. Also note that the pair-computations within "
		  "this module will use the cutoff distance defined by the used weighting "
		  "function."),
		(Node*) m_particleCache);

//	CPPUNIT_ASSERT_EQUAL
//  	(true, m_particleCache -> properties().attrExists("weightingFunction"));
//
//	CPPUNIT_ASSERT_EQUAL
//  	(string("Symbol name of the externally defined weighting function to be used. "
//  "SPH2ndDerivCorr will use the derivative of the defined weighting "
//  "(interpolation) function. Also note that the pair-computations within "
//  "this module will use the cutoff distance defined by the used weighting "
//  "function."),
//  		m_particleCache -> properties().attrByName("weightingFunction").description);

	CppunitHelper::testPropAttr(
		string("volume"),
		string(
				"Symbol name of the externally computed SPH particle volume."),
		(Node*) m_particleCache);

	CppunitHelper::testPropAttr(
		string("SPH1stDerivCorr"),
		string(
				"Symbol name of the externally computed correction matrix for the SPH-"
				    "discretisation of the 1st derivative."),
		(Node*) m_particleCache);

	CppunitHelper::testPropAttr(
		string("symbol"),
		string(
				"Symbol name of the correction matrix for the SPH-discretisation of the "
				  	"2nd derivative computed by this module."),
		(Node*) m_particleCache);

	CPPUNIT_ASSERT_EQUAL(string("B2"), m_particleCache -> mySymbolName());

}

void PCa2ndSPHDerivCorrTest :: testPropClassName()
{
	CPPUNIT_ASSERT_EQUAL
		(string("PCa2ndSPHDerivCorr"),
				m_particleCache -> returnProperties().className());
}

void PCa2ndSPHDerivCorrTest :: testSpeciesAttr()
{
	CppunitHelper::testPropAttr(
		string("species"),
		string(
				"Name for the species of the particles, this Symbol is used for. "
				  		"'species = \"ALL\"' is not allowed for this module."),
		(Node*) m_particleCache);
}

void PCa2ndSPHDerivCorrTest :: setupTest (void)
{
  // Currently we avoid testing for useful or reasonably faked xml-input.
	// Therefore we can do the same test as for the parent class
	// FIXME: find a way for a more in-depth test either by testing faked
	// xml-input or by faking the values of the variables associated to the xml-attributes
	ParticleCacheTest::setupTest();
}

// FIXME: tests missing:
// void PCa2ndSPHDerivCorrTest::initSystemMatricesTest (void)
// void PCa2ndSPHDerivCorrTest::buildSystemMatricesTest (void)
// void PCa2ndSPHDerivCorrTest::checkInputSymbolExistencesTest (void)
// if need be: void PCa2ndSPHDerivCorrTest::precomputeTest (void)

void PCa2ndSPHDerivCorrTest :: computeCacheForTest (void)
{
	CPPUNIT_ASSERT_EQUAL_MESSAGE(
			"Test can only run succesfully for code compiled with SPACE_DIMS = 3",
			3, SPACE_DIMS
	);

	size_t resultEntries = PCa2ndSPHDerivCorr::return_s_matEntries();

	const size_t colour = 0;

	FakePCa2ndSPHDerivCorr* fakePCa2ndSPHDerivCorr
		= (FakePCa2ndSPHDerivCorr*) m_particleCache;

	size_t offDiagIdx;

	// setup the DataFormat
	Particle::s_tag_format.resize(1);
	Particle::s_tag_format[colour].setFormatAndAlloc(new DataFormat());

	// Matrices are externally set by unittest, so no pair loop for building them
	fakePCa2ndSPHDerivCorr -> setPairLoopToDo(false);

	fakePCa2ndSPHDerivCorr -> setColour(colour);

	// Add required attributes

	size_t systemMatOffset = Particle::s_tag_format[colour].addAttribute(
			"testSysMat", DataFormat::VECTOR_DOUBLE, false, "testSysMat")
			.offset;
	fakePCa2ndSPHDerivCorr -> setSystemMatOffset(systemMatOffset);

	size_t resultOffset
		= Particle::s_tag_format[colour]
				.addAttribute("testResult", DataFormat::TENSOR, false, "testResult")
				.offset;
	fakePCa2ndSPHDerivCorr -> setOutputOffset(resultOffset);

	// set correct size for system matrix
  Particle::s_tag_format[colour].vectorDoubleByOffset(systemMatOffset)
  	-> resize
				(resultEntries * resultEntries);

  // create 1 Particle on which to solve the linear system
  Particle* p = new Particle(colour);

  vector<double>& sysMat
		= *((p->tag.vectorDoubleByOffset(systemMatOffset)).value());

  tensor_t& result = p->tag.tensorByOffset(resultOffset);

  ///// CASE 1: //////////////////

  // set values for system matrix

	// We just take the unit matrix such that solution = RHS
  for (size_t i = 0; i < resultEntries; ++i) {
    for (size_t j = 0; j < resultEntries; ++j) {
    	if(i == j)
    	 	sysMat[j + i * resultEntries] = 1.;
    	else sysMat[j + i * resultEntries] = 0.;
    }
  }

  fakePCa2ndSPHDerivCorr -> computeCacheFor(p);

  // check result
  offDiagIdx = SPACE_DIMS;
  for (size_t i = 0; i < SPACE_DIMS; ++i) {
    for (size_t j = i+1; j < SPACE_DIMS; ++j) {
    	// first argument is a RHS entry
    	CPPUNIT_ASSERT_EQUAL(
    		GSL_VECTOR_GET(
    			PCa2ndSPHDerivCorr::return_s_rhs(), i == j ? i : offDiagIdx
    						),
					result(i, j)
				);
    	++offDiagIdx;

    	// is the solution symmetric? Should be by coding
    	CPPUNIT_ASSERT_EQUAL(result(j, i), result(i, j));
    }
  }

  ///// CASE 2: //////////////////

  // set values for system matrix

  // Rouse matrix:
  // -2  1 0 ...
  //  1 -2 1 0 ...
  //  ...
  //         ... 0 1 -2
  for (size_t i = 0; i < resultEntries; ++i) {
    for (size_t j = 0; j < resultEntries; ++j) {
     	if(i == j)
     	 	sysMat[j + i * resultEntries] = -2.;
     	else if(i == j+1 || i == j-1)
     	 	sysMat[j + i * resultEntries] = 1.;
     	else
     	 	sysMat[j + i * resultEntries] = 0.;
    }
  }

  fakePCa2ndSPHDerivCorr -> computeCacheFor(p);

  // check result

  offDiagIdx = SPACE_DIMS;
  const double reference[6]
		= { 2.14286, 3.28571, 3.42857, 2.57143, 1.71429, 0.857143 };

  // preliminary test of test code
  assert(reference[3] == 2.57143);

    for (size_t i = 0; i < SPACE_DIMS; ++i) {

    	// check of diagonal entries (1st 3 in reference solution)
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(reference[i], result(i, i), 0.00001);

      for (size_t j = i+1; j < SPACE_DIMS; ++j) {

      	if(SPACE_DIMS == 3)
      		// the final off-diagonal matrix entry in result(i, j) is 1/2 of the
      		// result from the linear system
      		CPPUNIT_ASSERT_DOUBLES_EQUAL
						(0.5*reference[offDiagIdx], result(i, j), 0.00001);
      	++offDiagIdx;

      	// is the solution symmetric? Should be by coding
      	CPPUNIT_ASSERT_EQUAL(result(j, i), result(i, j));
      }
    }

    ///// CASE 3: NON-SYMMETRIC ///////

    // set values for system matrix

    // random non-symmetric (non-singular) system matrix
    const double
		nonSymMat[36]
			= {
				 1.52162, -2.51086, -1.80829, 2.38386, -1.55161, -1.98107,
				 0.427125, 1.17325, -0.725628, 1.24947, -1.15465, -0.0407065,
				 2.5236, -3.39753, 0.967568, 0.909269, -2.50996, -1.09852,
				 -0.0000895192, 0.977036, 0.22449, 0.81298, 0.787123, -0.788674,
				 -1.01879, 0.397601, -0.561835, -0.147591, -0.181743, -1.43793,
				 1.9561, 2.46243, 0.900747, -2.35028, -0.213386, 1.64877
				 };
    // set system matrix of particle
    for (size_t i = 0;
    		i < 36;
    		++i)
    	sysMat[i] = nonSymMat[i];

    // externally computed inverse
//    const double
//		nonSymMatInv[resultEntries * resultEntries]
//		 = {
//				 0.279972, -0.180977, -0.02525, 0.237885, -0.180358, 0.271603,
//				 -0.072233, 0.233535, -0.0225482, 0.18295, 0.121919, 0.0977924,
//				 -0.426695, 0.0879995, 0.341546, 0.336257, -0.0165198, -0.13652,
//				 -0.0761273, 0.296041, 0.03514, 0.260905, -0.32679, -0.220948,
//				 0.218461, -0.437125, -0.205361, 0.356229, -0.269867, 0.0499143,
//				 -0.0714143, 0.183279, -0.0994465, -0.321146, -0.459842, -0.0956825
//			 };

    // externally computed solution
    const double nonSymSol[6]
			= {
					-0.0737452, -0.138753, -0.00285126, -0.255053, 0.424026, -0.0124179
				 };

    fakePCa2ndSPHDerivCorr -> computeCacheFor(p);

    // check result
  	offDiagIdx = SPACE_DIMS;
    for (size_t i = 0; i < SPACE_DIMS; ++i) {
    	// check of diagonal entries (1st 3 in reference solution)
    	CPPUNIT_ASSERT_DOUBLES_EQUAL(nonSymSol[i], result(i, i), 0.00001);
      for (size_t j = i + 1; j < SPACE_DIMS; ++j) {
      	if(SPACE_DIMS == 3)
      		// the final off-diagonal matrix entry in result(i, j) is 1/2 of the
      		// result from the linear system
      		CPPUNIT_ASSERT_DOUBLES_EQUAL
						(0.5*nonSymSol[offDiagIdx], result(i, j), 0.00001);
      	++offDiagIdx;
      	// is the solution symmetric? Should be by coding
      	CPPUNIT_ASSERT_EQUAL(result(j, i), result(i, j));
      }
    }

}
