/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
 * David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
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



#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_IISPH_const_rho.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorIISPHconstRho> integrator_IISPH_const_rho("IntegratorIISPHconstRho");

const point_t IntegratorIISPHconstRho::dummyNullPoint = {0, 0, 0};



//---- Constructors/Destructor ----

IntegratorIISPHconstRho::IntegratorIISPHconstRho(Controller *controller):IntegratorPosition(controller)
{
  init();
}


IntegratorIISPHconstRho::~IntegratorIISPHconstRho()
{
}


//---- Methods ----

void IntegratorIISPHconstRho::init()
{
  // some modules need to know whether there is an Integrator,
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorIISPHconstRho");
  
  m_properties.setDescription
    ("Implicit incompressible SPH Integrator based on the algorithm proposed by Ihmsen at al. (IEEE Transactions on Visualization and Computer Graphics 20, 426 (2014)).\n"
     "NOTE THE FOLLOWING PROPERTIES, ASSUMPTIONS, LIMITATIONS, AND REQUIREMENTS:\n"
     "- If you use more than one species, read the description of attribute 'advVelocityName' below carefully." 
     "- This Integrator expects a Symbol storing a (generalised) mass for EVERY species in the simulation, no matter if the incompressible species actually interacts with a given species. In the simplest case, this Symbol just stores a constant mass for every particle of the given species, but it must exist."
     "- This Integrator enforces incompressibility within a SINGLE species. Multi-phase incompressibility is not yet (2017-01-18) supported.\n"
     "- Hard collisions with walls or similar objects are currently (2017-01-05) not supported. Using this integrator, such collisons will not be detected and particles may penetrate walls leading to loss of particles or to their motion in 'forbidden' areas.\n"
     "- All particles j from other species than the incompressible one (with particles i) excert the pressure Pi (!) on particle i. I.e. the particles j do not possess their own separately computed pressure Pj! For particles j forming rigid objects this should work. For particles j representing soft-matter, we should use a (on 2017-01-10 not yet implemented) multi-phase incompressible Integrator.\n"
     "- This Integrator does generally not compute pressure values for other species than the incompressible one.\n"
     "The velocity of all other species that the incompressible species is interacting with must be computed outside of IntegratorIISPHconstRho. For walls it might be zero (or s.th. that better fulfills e.g. no-slip if needed). For a moving rigid-body there will be a rigid-body Integrator and one still (2017-01-11) has to think how to do the coupling consistently"
     );

  STRINGPC
    (weightingFunction, m_kernelName,
     "Defines the weighting function (interpolation kernel) to be used. Note that for "
     "each pair of species, "
     "the maximum of the cutoff defined in the weighting function specified here and the cutoff used "
     "elsewhere in the simulation setup is applied.");     
  m_kernelName = "default";

  STRINGPC
    (otherSpecies, m_speciesList,
     "All the other species for which the interaction with the incompressible species (defined in attribute 'species' should be taken into account in the incompressibility algorithm. Multiple species should be separated by the pipe ('|') symbol. Additional Integrators introducing additional degrees of freedom for the incompressible species or any of the other interacting species must be defined BELOW IntegratorIISPHconstRho.");
  m_speciesList = "---";

  STRINGPC
    (genMassName, m_genMassName,
     "Name of the particle attribute storing the generalised particle masses of all species"
     );
  m_genMassName = "UNDEFINED";
  
  STRINGPC
    (advVelocityName, m_vAdvName,
     "Name of the attribute storing the velocity advected by the non-pressure forces.\n"
     "NOTE: This Integrator will compute the advected velocity for the incompressible species. "
     "For all other species, this Integrator will initialise the advected velocity once with 0. "
     "This means that the advected velocity for all other species will stay at this value "
     "forever unless the user introduces other modules to compute a Symbol with this name (in "
     "overwrite-mode). The latter might be necessary depending on the boundary conditions that "
     "should be implemented.");
  m_vAdvName = "IISPHvAdv";

  STRINGPC
    (pressureName, m_pressureName,
     "Name of the particle attribute storing the final value computed pressure."
     );
  m_pressureName = "UNDEFINED";

  STRINGPC
    (densityName, m_densityName,
     "Name of the particle attribute storing the density variable obeying incompressibility."
     );
  m_densityName = "UNDEFINED";
  
  DOUBLEPC
    (maxDensityError, m_maxDensityError, 0.,
     "Maximal maximum relative density error for determination of the fulfillment of the incompressibility "
     "condition. The maximum relative density error for a particle i is defined by \n"
     "max_i (abs(rho_i - rho0)/rho0 ) \n"
     "where rho_i is the local particle density and rho0 is the prescribed global reference density."
     );
  m_maxDensityError = 0.01;
  
  DOUBLEPC
    (avgDensityError, m_avgDensityError, 0.,
     "Maximal average relative density error for determination of the fulfillment of the incompressibility "
     "condition. The average relative density error for a particle i is defined by \n"
     "(1/N)*sum_i abs(rho_i - rho0)/rho0 \n"
     "where N is the total number of particles which should fulfil incompressibility, rho_i is the local "
     "particle density and rho0 is the prescribed global reference density."
     );
  m_avgDensityError = 0.005;

  DOUBLEPC
    (omega, m_omega, 0.,
     "Relaxation parameter for the relaxed Jacobi iteration\n"
     "Pnew_i = (1-omega)*Pold_i + omega*(RHS_i-a_ij*Pold_j)/a_ii\n"
     "where Pnew_i and Pold_i are the new and old pressure iteration values for particle i. The remaining "
     "parameters refer to the system of equations a_ij*P_j=RHS_i (summation over repeated indices) which "
     "is in fact solved. Admissible range is 0 < omega < 1."
     );
  m_omega = 0.5;
  
  DOUBLEPC
    (rho0, m_rho0, 0.,
     "Reference density for the incompressibility condition."
     );
  m_rho0 = HUGE_VAL;
  
  INTPC
    (nIterSteps, m_lMax, 0,
     "Maximum number of iteration steps for pressure computation. If it is exceeded, the program is "
     "aborted with an error message."
     );
  m_lMax = 10;
  
}

void IntegratorIISPHconstRho::setup()
{
  
  IntegratorPosition::setup();

  DataFormat::attribute_t tempAttr;

  if(m_rho0 == HUGE_VAL)
    throw gError("IntegratorIISPHconstRho::setup", "Attribute 'rho0' was not defined!");

  if(m_rho0 <= 0.)
    throw gError("IntegratorIISPHconstRho::setup", "Invalid value \"" + ObjToString(m_rho0) + "\" for attribute 'rho0'. Must be >0!");

  if(m_omega <= 0.)
    throw gError("IntegratorIISPHconstRho::setup", "Invalid value \"" + ObjToString(m_omega) + "\" for attribute 'omega'. Must be 0 < omega < 1!");

  if(m_omega >= 1.)
    throw gError("IntegratorIISPHconstRho::setup", "Invalid value \"" + ObjToString(m_omega) + "\" for attribute 'omega'. Must be 0 < omega < 1!");

  if(m_maxDensityError <= 0.)
    throw gError("IntegratorIISPHconstRho::setup", "Invalid value \"" + ObjToString(m_maxDensityError) + "\" for attribute 'maxDensityError'. Must be >0!");
  
  if(m_avgDensityError <= 0.)
    throw gError("IntegratorIISPHconstRho::setup", "Invalid value \"" + ObjToString(m_avgDensityError) + "\" for attribute 'avgDensityError'. Must be >0!");
  
  if(m_lMax < 1)
    throw gError("IntegratorIISPHconstRho::setup", "Invalid value \"" + ObjToString(m_lMax) + "\" for attribute 'nIterSteps'. Must be >0!");
  
  // interpolation kernel
  m_kernel = M_SIMULATION->findWeightingFunction(m_kernelName);
  
  // This simplifies some loops in the algorithm
  if(m_colour != 0)
    throw gError("IntegratorIISPHconstRho::setup", "The incompressible species was not defined as the very first species but this is required by IntegratorIISPHconstRho. To achieve this, no other modules (such as Integrators) should define a species before IntegratorIISPHconstRho.");


  // -------- START: Create those other colours we interact with ----------------------
  
  if(m_speciesList == "---")
    MSG_DEBUG("IntegratorIISPHconstRho::setup", "User selection for other species: \"---\". So no other species will be created and considered.");
  else {
    string remainingSpeciesList = m_speciesList;
    bool run = true;
    size_t colourTester = 0;
    
    while (run) {
      string cur;
      size_t pos = remainingSpeciesList.find('|');
      
      if (pos == string::npos) {
	run = false;
	cur = remainingSpeciesList;
      } 
      else {
	cur = string(remainingSpeciesList, 0, pos);
	remainingSpeciesList = string(remainingSpeciesList, pos+1);
      }
      
      ++colourTester;
      // adding colour and verifying that it has the right number
      assert(colourTester == Integrator::getColourAndAdd(cur));
      MSG_DEBUG("IntegratorIISPHconstRho::setup", "Added new species " << cur << " with colour index " << colourTester << ".");
      
    } // end of while (run)
  } // end of else of if(m_speciesList == "---")
  
  // -------- END: Create those other colours we interact with ----------------------
  

  // set cutoffs, activate pairs, and
  // create the list of ColourPairs used here and containing m_colour ONCE AND ONLY ONCE
  // at this point only this Integrator has created colours, hence fine to loop over all ColourPairs
  FOR_EACH_COLOUR_PAIR
    (M_MANAGER,
     // Should not happen if m_colour == 0, right?!?
     if(cp->firstColour() != m_colour) assert(cp->secondColour() != m_colour);

     // set cutoff and activate pairs
     // might very well be that the cutoff is increased later by other modules for some ColourPairs
     cp->setCutoff(m_kernel->cutoff());
     cp->setNeedPairs(true);
     
     // list mixed ColourPairs
     if(cp->firstColour() == m_colour && cp->secondColour() != m_colour) {
       m_mixedColourPairs.push_back(cp);
       MSG_DEBUG("IntegratorIISPHconstRho::setup", "I will take into account mixed ColourPair (" << cp->firstColour() << " ," << cp->secondColour() << ").");
     }

     );

  
  // --- START: internal Symbol for advected velocity (not yet existing and to be created) ------------------
  
  // Advected velocity Symbol for this colour  
  assert(m_vAdvOffset.size() == 0);
  // m_vAdvOffset.push_back(addNewSymbol(m_colour, "__IISPHvAdv"));
  m_vAdvOffset.push_back(addNewAttr(m_colour, m_vAdvName, DataFormat::POINT));

  // Create advected velocity Symbol for all second colours of the mixed ColourPairs
  for(vector<ColourPair*>::iterator cpit = m_mixedColourPairs.begin(); cpit != m_mixedColourPairs.end(); ++cpit)
    {
      ColourPair* cp = (*cpit);
      assert(m_colour == cp -> firstColour());
      size_t secondColour = cp -> secondColour();

      m_vAdvOffset.push_back(addNewAttr(secondColour, m_vAdvName, DataFormat::POINT));
      
    }

  // --- END: internal Symbol for advected velocity (not yet existing and to be created) --------------------

  
  // --- START: generalised mass Symbol PART 1 ----------------------------------

  if(m_genMassName == "UNDEFINED")
    throw gError("IntegratorIISPHconstRho::setup", "Attribute genMassName has value \"UNDEFINED\".");
  // reserve memory; at this stage only this Integrator has created colours, so the call to Phase is fine
  m_genMassOffset.resize(M_PHASE->nColours());
  
  // we can not assign offsets here because the mass attributes are created externally and after this setup
  // see IntegratorIISPHconstRho::isAboutToStart()
  
  // --- END: generalised mass Symbol PART 1 ----------------------------------

  
  // add internal Symbols and get offset; hence Symbols should not yet exist
  
  m_diiOffset = addNewAttr(m_colour, "__IISPHdii", DataFormat::POINT);
  m_aiiOffset = addNewAttr(m_colour, "__IISPHaii", DataFormat::DOUBLE);
  m_advDensityOffset = addNewAttr(m_colour, "__IISPHadvDensity", DataFormat::DOUBLE);
  m_pforcePairIncrOffset = addNewAttr(m_colour, "__IISPHpIncr", DataFormat::POINT);
  m_pressureIterOffsetOld = addNewAttr(m_colour, "__IISPHpIterSlot1", DataFormat::DOUBLE);
  m_pressureIterOffsetNew = addNewAttr(m_colour, "__IISPHpIterSlot2", DataFormat::DOUBLE);
  m_precomputeOffset = addNewAttr(m_colour, "__IISPHprecompute", DataFormat::DOUBLE);
  m_pforceIncrOffset = addNewAttr(m_colour, "__IISPHpForceIncr", DataFormat::POINT);
  m_iterDensityOffset = addNewAttr(m_colour, "__IISPHiterDensity", DataFormat::DOUBLE);

  // pressure to be added with user defined name
  if(m_pressureName == "UNDEFINED")
    throw gError("IntegratorIISPHconstRho::setup", "Attribute pressureName has value \"UNDEFINED\".");
  m_pressureOffset = addNewAttr(m_colour, m_pressureName, DataFormat::DOUBLE);
  
  // the density Symbol
  if(m_densityName == "UNDEFINED")
    throw gError("IntegratorIISPHconstRho::setup", "Attribute densityName has value \"UNDEFINED\".");  
  m_densityOffset = addNewAttr(m_colour, m_densityName, DataFormat::DOUBLE);
  
}


// FIXME: looks like a method that should go, e.g., into class DataFormat or Data
size_t IntegratorIISPHconstRho::addNewAttr(size_t colour, string symbolName, DataFormat::datatype_t datatype, bool persistency/* = true*/)
{
  if(Particle::s_tag_format[colour].attrExists(symbolName))
    throw gError("IntegratorIISPH::addNewSymbol", ": Symbol " + symbolName + " is already existing for species '" + M_MANAGER->species(colour) + "'. If you have chosen this name for one of your other Symbols, choose a different one! If you haven't then this is an internal error and you should file a bug report.");

  // OK, we can create it and return the offset
  return Particle::s_tag_format[colour].addAttribute(symbolName, datatype, persistency, symbolName).offset;
}


// FIXME: Can this be merged with setupAfterParticleCreation() ?
void IntegratorIISPHconstRho::isAboutToStart()
{
  Phase *phase = M_PHASE;
  m_dt = M_CONTROLLER->dt();
  DataFormat::attribute_t tempAttr;
  
  // initialise ALL advected velocities to zero. See help-text about the attribute name
  for (size_t colour = 0; colour < phase->nColours(); ++colour) { 
    FOR_EACH_PARTICLE_C
      (phase, colour,
       for(size_t i = 0; i < SPACE_DIMS; ++i)
	 __iSLFE->tag.pointByOffset(m_vAdvOffset[colour])[i] = 0.;
       );

  }


  // --- START: generalised mass Symbol PART 2 ----------------------------------
  // see IntegratorIISPHconstRho::setup for PART 1
  
  // assign offset of generalised mass for m_colour; other colours follow below
  tempAttr = Particle::s_tag_format[m_colour].attrByName(m_genMassName);

  // consistency check
  if(tempAttr.datatype != DataFormat::DOUBLE)
    throw gError("IntegratorIISPHconstRho::isAboutToStart", "the symbol " + m_genMassName +
		 " is registerd as a non-scalar for species " +
		 M_MANAGER->species(m_colour));

  // assign offsets
  m_genMassOffset[m_colour] = tempAttr.offset;
  
  // assign offset of generalised mass for all second colours of mixed ColourPairs
  for(vector<ColourPair*>::iterator cpit = m_mixedColourPairs.begin(); cpit != m_mixedColourPairs.end(); ++cpit)
    {
      ColourPair* cp = (*cpit);
      assert(m_colour == cp -> firstColour());
      size_t secondColour = cp -> secondColour();
      
      tempAttr = Particle::s_tag_format[secondColour].attrByName(m_genMassName);
      
      // consistency check
      if(tempAttr.datatype != DataFormat::DOUBLE)
	throw gError("IntegratorIISPHconstRho::isAboutToStart", "the symbol " + m_genMassName +
		     " is registerd as a non-scalar for species " +
		     cp->manager()->species(secondColour));

      // assign offsets
      m_genMassOffset[secondColour] = tempAttr.offset;

    }

  // --- END: generalised mass Symbol PART 2 ----------------------------------

  // initialise forces to zero  
  size_t counter = 0;
  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     for (int j = 0; j < FORCE_HIST_SIZE; j++)
       __iSLFE->force[j].assign(0);
     ++counter;
    );
  if(counter == 0)
    throw gError("IntegratorIISPHconstRho::isAboutToStart", "no free particles found for species " + m_species + "! Don't instantiate an Integrator for positions and velocities in that case. Use another module to create the species.");
  // FIXME: so we need some SpeciesCreator to make it more transparent
  // FIXME: put all in this function into the general setup for Nodes after the particle creation or into s.th. even more general

}


void IntegratorIISPHconstRho::integrateStep1()
{
  M_PHASE->invalidatePositions((IntegratorPosition*) this);
}


void IntegratorIISPHconstRho::integrateStep2()
{
  ColourPair* cpSelf = M_MANAGER->cp(m_colour, m_colour);
  Phase *phase = M_PHASE;
  size_t force_index = M_CONTROLLER->forceIndex();
  size_t genMassOffsetSelf = m_genMassOffset[m_colour];
  size_t vAdvOffsetSelf = m_vAdvOffset[m_colour];
  double dtsq = m_dt*m_dt;
  // const point_t dummyNullPoint = {0, 0, 0};
  
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handle
     const Data& pTag = i->tag;
     
     // initialisation for below
     pTag.doubleByOffset(m_aiiOffset) = 0.;
     for(size_t j = 0; j < SPACE_DIMS; ++j)
       pTag.pointByOffset(m_diiOffset)[j] = 0.;
     
     // the initial density disturbed below by the non-pressure forces, is the previous density,
     // but here we initialise with zero and add the previous density below.
     pTag.doubleByOffset(m_advDensityOffset) = 0.;
       
     // advected velocity
     pTag.pointByOffset(vAdvOffsetSelf) = i->v + m_dt*i->force[force_index]/(pTag.doubleByOffset(genMassOffsetSelf));

     );
  
  // --- START: diagonal part of d-matrix ----------------------------------------

  // NOTE: The summation over further coulours/species is only performed for the precomputed matrix-elements dii and not for the on-the-fly computed elements dij. This implies the ASSUMPTION that all particles j from other species than the incompressible one excert the pressure Pi (!) on particle i. I.e. the particles j do not possess their own separately computed pressure Pj! For particles j forming rigid objects this should work. For particles j representing soft-matter, we should use a (on 2017-01-10 not yet implemented) multi-phase incompressible Integrator. 
  
  // the self pair-specific part
  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorIISPHconstRho, cpSelf,
     // handles
     const Data& p1Tag = pair->firstPart()->tag;
     const Data& p2Tag = pair->secondPart()->tag;
     
     point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);

     // first particle
     p1Tag.pointByOffset(m_diiOffset) +=
     // this may just be the mass or a corrected effective mass;
     // the user must decide by precomputation what to take for each species
     p2Tag.doubleByOffset(genMassOffsetSelf)
     // the gradient of the interpolation function
     *interpolGradient;

     // "-=" for 2nd particle since antisymmetric
     p2Tag.pointByOffset(m_diiOffset) -=
     // this may just be the mass or a corrected effective mass;
     // the user must decide by precomputation what to take for each species
     p1Tag.doubleByOffset(genMassOffsetSelf)
     // the gradient of the interpolation function
     *interpolGradient;

     );   
  
  // the non-self pair-specific part
    for(vector<ColourPair*>::iterator cpit = m_mixedColourPairs.begin(); cpit != m_mixedColourPairs.end(); ++cpit)
    {
      ColourPair* cp = (*cpit);
      size_t genMassOffset = m_genMassOffset[cp->secondColour()];
      FOR_EACH_PAIR__PARALLEL
	(IntegratorIISPHconstRho, cp,
	 // setup should have checked that fluid species is the very first,
	 // so the following call is the only possible one
	 pair->firstPart()->tag.pointByOffset(m_diiOffset) +=
	 // this may just be the mass or a corrected effective mass;
	 // the user must decide by precomputation what to take for each species
	 pair->secondPart()->tag.doubleByOffset(genMassOffset)
	 // the gradient of the interpolation function
	 *(-pair->cartesian()*m_kernel->weight(pair, dummyNullPoint));
	 );
    }

  // the particle-specific prefactor
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // density must have been computed before, either externally or
     // in the previous time step further down in this method
     const double& pDensity = i->tag.doubleByOffset(m_densityOffset);
     i->tag.pointByOffset(m_diiOffset) *= -dtsq/(pDensity*pDensity);
     );
    
  // --- END: diagonal part of d-matrix ----------------------------------------

  // --- START: diagonal part of a-matrix ----------------------------------------

  // No need for other ColourPairs due to currently made ASSUMPTIONS:
  // Basically we compute matrix entries for incompressibility WITHIN a SINGLE species  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorIISPHconstRho, cpSelf,
     // handles
     const Data& p1Tag = (pair->firstPart()->tag);
     const Data& p2Tag = (pair->secondPart()->tag);
     const double& p1Mass = p1Tag.doubleByOffset(genMassOffsetSelf);
     const double& p2Mass = p2Tag.doubleByOffset(genMassOffsetSelf);
     const double& p1Density = p1Tag.doubleByOffset(m_densityOffset);
     const double& p2Density = p2Tag.doubleByOffset(m_densityOffset);
     // the gradient of the interpolation function
     point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);

     p1Tag.doubleByOffset(m_aiiOffset) +=
     p2Mass*(
	     p1Tag.pointByOffset(m_diiOffset)
	     // second term is -dji
	     // minus-signs: 1. -dtsq^2, 2. gradWji=-gradWij, 3. overall "-" in sum 
	     - dtsq*p1Mass*interpolGradient/(p1Density*p1Density)
	     )*interpolGradient; // should be a scalar product
     
     p2Tag.doubleByOffset(m_aiiOffset) -=
     p1Mass*(
	     p2Tag.pointByOffset(m_diiOffset)
	     + dtsq*p2Mass*interpolGradient/(p2Density*p2Density)
	     )*interpolGradient; // should be a scalar product; minus sign is in the "-=" above

     // --- END: diagonal part of a-matrix ----------------------------------------

     
     // --- START: advected density due to non-pressure forces - self-pair-part --------

     // within the pair loop we only compute the increment (without factor m_dt)
     double deltaVelAdvDotInterpolGradient =
     (p1Tag.pointByOffset(vAdvOffsetSelf) - p2Tag.pointByOffset(vAdvOffsetSelf))*interpolGradient;

     // two "+=" because symmetric expression
     p1Tag.doubleByOffset(m_advDensityOffset) += p2Mass*deltaVelAdvDotInterpolGradient;
     p2Tag.doubleByOffset(m_advDensityOffset) += p1Mass*deltaVelAdvDotInterpolGradient;

     MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "p1Mass(" << pair->firstPart()->mySlot << ")=" << p1Mass << ", p2Mass(" << pair->secondPart()->mySlot << ")=" << p2Mass<< ", p1vAdv=" << p1Tag.pointByOffset(vAdvOffsetSelf) << ", p2vAdv=" << p2Tag.pointByOffset(vAdvOffsetSelf) << ", interpolGradient=" << interpolGradient << ", pair->r=" << pair->cartesian() << ", -W'(r)/r=" << m_kernel->weight(pair, dummyNullPoint) << ", eij=" << pair->cartesian()/pair->abs());
	  
     // if(pair->firstPart()->mySlot==48 || pair->secondPart()->mySlot==48)
        MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "advDensityPairLoop(" << pair->firstPart()->mySlot << "," << pair->secondPart()->mySlot << "): incr1=" << p2Mass*deltaVelAdvDotInterpolGradient << ", incr2=" << p1Mass*deltaVelAdvDotInterpolGradient << ", rhoAdvNew1=" << p1Tag.doubleByOffset(m_advDensityOffset) << ", rhoAdvNew2=" << p2Tag.doubleByOffset(m_advDensityOffset));
     
     ); // end FOR_EACH_PAIR__PARALLEL

  // --- END: advected density due to non-pressure forces - self-pair-part --------------


  // --- START: advected density due to non-pressure forces - non-self-pair-part --------

    for(vector<ColourPair*>::iterator cpit = m_mixedColourPairs.begin(); cpit != m_mixedColourPairs.end(); ++cpit)
    {
      ColourPair* cp = (*cpit);
      size_t secondColour = cp->secondColour();
      size_t genMassOffset2ndColour = m_genMassOffset[secondColour];
      FOR_EACH_PAIR__PARALLEL
	(IntegratorIISPHconstRho, cp,
	 // handles
	 const Data& p1Tag = (pair->firstPart()->tag);
	 const Data& p2Tag = (pair->secondPart()->tag);
	 const double& p2Mass = p2Tag.doubleByOffset(genMassOffset2ndColour);

	 point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);

	 // within the pair loop we only compute the increment (without factor m_dt)
	 // ASSUMPTION: the second velocity must be computed outside of IntegratorIISPHconstRho;
	 // for walls it might be zero (or s.th. that better fulfills e.g. no-slip if needed);
	 // for a moving rigid-body there will be a rigid-body Integrator and one still (2017-01-11)
	 // has to think how to do the coupling consistently
	 double deltaVelAdvDotInterpolGradient =
	 (p1Tag.pointByOffset(m_vAdvOffset[m_colour]) - p2Tag.pointByOffset(m_vAdvOffset[secondColour]))
	 *interpolGradient;

	 // two "+=" because symmetric expression
	 p1Tag.doubleByOffset(m_advDensityOffset) += p2Mass*deltaVelAdvDotInterpolGradient;
	 // no advected density for 2nd species needed so far (2017-01-11)
	 // p2Tag.doubleByOffset(m_advDensityOffset) += p1Mass*deltaVelAdvDotInterpolGradient;
	 
	 );
    }

  // --- END: advected density due to non-pressure forces - non-self-pair-part --------

  
  // --- START: advected density due to non-pressure forces - particle-wrapup-part --------
  
  // the initial value for the summation is the previous density, which is added now
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // but first the missing factor m_dt
     i->tag.doubleByOffset(m_advDensityOffset) *= m_dt;

     // if(i->mySlot == 48) {
       // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "pSlot=" << i->mySlot << ", rhoAdvBeforeFinal(=DifferenceToOldRho)=" << i->tag.doubleByOffset(m_advDensityOffset) );
       // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "pSlot=" << i->mySlot << ", OldRho-rho0=" << i->tag.doubleByOffset(m_densityOffset)-m_rho0 );
     // }
     
     i->tag.doubleByOffset(m_advDensityOffset) += i->tag.doubleByOffset(m_densityOffset);

     // EXTENSION: here comes the line P_i^m=0 = P_i(t-dt) in IntegratorIISPH_EOS

     // if(i->mySlot == 48)
       MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "pSlot=" << i->mySlot << ", final rhoAdv=" << i->tag.doubleByOffset(m_advDensityOffset) << ", ...-rho0=" << i->tag.doubleByOffset(m_advDensityOffset)-m_rho0);
     
     );
  
  // --- END: advected density due to non-pressure forces - particle-wrapup-part --------

  
  // --- START: Relaxed Jacobi iteration ----------------------------------------

  // initialisation of iteration
  size_t l = 0;
  double deltaDensityMax = HUGE_VAL;
  double deltaDensityAvg = HUGE_VAL;

  // initial value for pressure iteration
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handles
     const Data& pTag = i->tag;
     // factor of 0.5 is best initial value according to Ihmsen et al.
     // (IEEE Transactions on Visualization and Computer Graphics 20, 426 (2014))
     pTag.doubleByOffset(m_pressureIterOffsetNew) = 0.5*pTag.doubleByOffset(m_pressureOffset);
     );

  
  // --- START: iterated pair-pressure increment D_i ----------------------------------------

  // we compute initial value here and updated values in while loop after having new pressures
  computePairPressureIncrement();            

  // --- END: iterated pair-pressure increment D_i ----------------------------------------

  // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "BEFORE ITER: deltaDensityMax = " << deltaDensityMax << ", m_maxDensityError = " << m_maxDensityError << ", deltaDensityAvg = " << deltaDensityAvg << ", m_avgDensityError = " << m_avgDensityError);
  
  // The iteration
  while(deltaDensityMax > m_maxDensityError || deltaDensityAvg > m_avgDensityError || l < 2)
    {
      // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "ITER-START, l = " << l);

      // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "IN ITER: deltaDensityMax = " << deltaDensityMax << ", m_maxDensityError = " << m_maxDensityError << ", deltaDensityAvg = " << deltaDensityAvg << ", m_avgDensityError = " << m_avgDensityError);
      
      if(l > m_lMax)
	throw gError("IntegratorIISPHconstRho::integrateStep2", "Iteration steps exceeded user defined maximum of " + ObjToString(m_lMax) + " without convergence. Here are the error limits and errors:\n"
		     "maxDensityErrorLimit = " + ObjToString(m_maxDensityError) +
		     " (current = " + ObjToString(deltaDensityMax) + ")\n"
		     "avgDensityErrorLimit = " + ObjToString(m_avgDensityError) +
		     " (current = " + ObjToString(deltaDensityAvg) + ")\n"
		     "\nAborting.");
      
      deltaDensityMax = 0.;
      deltaDensityAvg = 0.;
      
      // --- START: new pressure value ----------------------------------------

      // Set newest pressure to old; set previous old pressure to new for new summation
      size_t temp = m_pressureIterOffsetOld;
      m_pressureIterOffsetOld = m_pressureIterOffsetNew;
      m_pressureIterOffsetNew = temp;

      // reset to zero for new summation (we will be working on the pressureIterOffsetNew slot!!!)
      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(phase, m_colour, this,
	 // this will be used for the next pressure iteration and is hence reset here
	 i->tag.doubleByOffset(m_pressureIterOffsetNew) = 0.;
	 );
      
      // self-pair contribution
      
      FOR_EACH_PAIR__PARALLEL
	(IntegratorIISPHconstRho, cpSelf,
	 
	 // handles
	 const Data& p1Tag = (pair->firstPart()->tag);
	 const Data& p2Tag = (pair->secondPart()->tag);
	 const double& p1Mass = p1Tag.doubleByOffset(genMassOffsetSelf);
	 const double& p2Mass = p2Tag.doubleByOffset(genMassOffsetSelf);
	 const point_t& p1Dii = p1Tag.pointByOffset(m_diiOffset);
	 const point_t& p2Dii = p2Tag.pointByOffset(m_diiOffset);
	 const point_t& p1PairPincr = p1Tag.pointByOffset(m_pforcePairIncrOffset);
	 const point_t& p2PairPincr = p2Tag.pointByOffset(m_pforcePairIncrOffset);
	 const double& p1PressureIterOld = p1Tag.doubleByOffset(m_pressureIterOffsetOld);
	 const double& p2PressureIterOld = p2Tag.doubleByOffset(m_pressureIterOffsetOld);
	 // previously computed in computePairPressureIncrement()
	 const double& p1PrecomputedNegDtSq_M_Piter_divRhoSq = p1Tag.doubleByOffset(m_precomputeOffset);
	 const double& p2PrecomputedNegDtSq_M_Piter_divRhoSq = p2Tag.doubleByOffset(m_precomputeOffset);
	 point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);
	 
	 // We compute -Sum_j m_j * [ D_i - d_jj*P_j_old - D_j + d_ji*P_i_old ] * nablaW_ij
	 // with
	 // d_ji*P_i_old = -(dt^2*m_i*P_i_old/rho_i^2)*nablaW_ji
	 // = -(dt^2*m_i*P_i_old/rho_i^2)*(-nablaW_ij)
	 // = (-1)(-(dt^2*m_i*P_i_old/rho_i^2))*(nablaW_ij)
	 // = -p1PrecomputedNegDtSq_M_Piter_divRhoSq*interpolGradient
	 // and
	 // D_i = p1Tag.pointByOffset(m_pforcePairIncrOffset) = p1PairPincr

	  // if(pair->firstPart()->mySlot == 48 || pair->secondPart()->mySlot == 48)
	  //   {
	  //     MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "newPpairLoop: Pbefore(" << pair->firstPart()->mySlot << ")=" << p1Tag.doubleByOffset(m_pressureIterOffsetNew));
	  //     MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "newPpairLoop: Pbefore(" << pair->secondPart()->mySlot << ")=" << p2Tag.doubleByOffset(m_pressureIterOffsetNew));
	  //   }
	     
	 p1Tag.doubleByOffset(m_pressureIterOffsetNew) -=
	 p2Mass*(
		 p1PairPincr - p2Dii*p2PressureIterOld
		 - p2PairPincr - p1PrecomputedNegDtSq_M_Piter_divRhoSq*interpolGradient
		 )*interpolGradient;

	 // last term in the bracket is
	 // d_ij*P_j_old = -(dt^2*m_j*P_j_old/rho_j^2)*nablaW_ij
	 // = p2PrecomputedNegDtSq_M_Piter_divRhoSq*interpolGradient
	 p2Tag.doubleByOffset(m_pressureIterOffsetNew) +=
	 p1Mass*(
		 p2PairPincr - p1Dii*p1PressureIterOld
		 - p1PairPincr + p2PrecomputedNegDtSq_M_Piter_divRhoSq*interpolGradient
		 )*interpolGradient;

	  // if(pair->firstPart()->mySlot == 48 || pair->secondPart()->mySlot == 48)
	  //   {
	  //     MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "newPpairLoop: Pafter(" << pair->firstPart()->mySlot << ")=" << p1Tag.doubleByOffset(m_pressureIterOffsetNew));
	  //     MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "newPpairLoop: Pafter(" << pair->secondPart()->mySlot << ")=" << p2Tag.doubleByOffset(m_pressureIterOffsetNew));
	  //     MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "newPpairLoop(" << pair->firstPart()->mySlot << "," << pair->secondPart()->mySlot << "): p1Mass = " << p1Mass << ", p2Mass = " << p2Mass << ", p1PairPincr = " << p1PairPincr << ", p2PairPincr = " << p2PairPincr << ", p1Dii = " << p1Dii << ", p2Dii = " << p2Dii << ", p1PressureIterOld = " << p1PressureIterOld << ", p2PressureIterOld = " << p2PressureIterOld << ", p1PrecomputedNegDtSq_M_Piter_divRhoSq = " << p1PrecomputedNegDtSq_M_Piter_divRhoSq << ", p2PrecomputedNegDtSq_M_Piter_divRhoSq = " << p2PrecomputedNegDtSq_M_Piter_divRhoSq << ", interpolGradient = " << interpolGradient);
	  //   }
	 
	 );

      // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "newPselfPairLoop FINISHED");
      
      // non-self-pair contribution
      
      for(vector<ColourPair*>::iterator cpit = m_mixedColourPairs.begin(); cpit != m_mixedColourPairs.end(); ++cpit)
	{
	  ColourPair* cp = (*cpit);
	  size_t genMassOffsetSecondColour = m_genMassOffset[cp->secondColour()];
	  FOR_EACH_PAIR__PARALLEL
	    (IntegratorIISPHconstRho, cp,
	     // handles
	     const Data& p1Tag = (pair->firstPart()->tag);
	     const Data& p2Tag = (pair->secondPart()->tag);
	     const point_t& p1PairPincr = p1Tag.pointByOffset(m_pforcePairIncrOffset);
	     const double& p2Mass = p2Tag.doubleByOffset(genMassOffsetSecondColour);
	     // the gradient of the interpolation function
	     point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);

	     // we compute -sum_b [m_b * nablaW_ib . sum_j [d_ij * P_j_old] ]
	     // =          -sum_b [m_b * nablaW_ib . D_i ]
	     // where b is a particle from the 2nd species
	     // and D_i = p1PairPincr
	     // EXTENSION: m_b might become m_b(rho_eq_i) in IntegratorIISPHEOS
	     p1Tag.doubleByOffset(m_pressureIterOffsetNew) -=
	     // last product is scalar product
	     p2Mass*interpolGradient*p1PairPincr;

	     // ASSUMPTION: no need for an update of p2 yet (2017-01-12)
	     
	     );
	}

      // particle wrap-up-contribution

      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(phase, m_colour, this,
	 // handles
	 const Data& pTag = i->tag;
	 double& pressureIterNew = pTag.doubleByOffset(m_pressureIterOffsetNew);
	 const double& pressureIterOld = pTag.doubleByOffset(m_pressureIterOffsetOld);
	 const double& rhoAdv = pTag.doubleByOffset(m_advDensityOffset);
	 const double& aii = pTag.doubleByOffset(m_aiiOffset);

	 // if(i->mySlot == 48)
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "New Piter: PnewTemp = " << pressureIterNew << ", m_omega = " << m_omega << ", aii=" << aii << ", pIndex=" << i->mySlot);

	 
	 // FIXME: the difference m_rho0 - rhoAdv will initially often be that small
	 // that it numerically ends up to be zero, since we compute the difference of
	 // two "big numbers". Think if we shouldn't compute this
	 // difference in such a way that the small difference survives, e.g., by
	 // computing(which we do anyway) and directly using the increment
	 // deltaRhoAdv = rhoAdv - rhoOld. We would then also need the accumulated
	 // increment accDeltaRho = rhoNewFinal - rho0 = deltaRhoThisStep + rhoOld - rho0.
	 // Summing the two equations gives (DOUBLE-CHECK THIS CALCULATION!)
	 // rho0 - rhoAdv = deltaRhoThisStep - deltaRhoAdv - accDeltaRho.
	 // So the RHS is computing exclusively with small deltas ! 
	 pressureIterNew = (1-m_omega)*pressureIterOld
	 + m_omega*(m_rho0 - rhoAdv + pressureIterNew)/aii;

	 // if(i->mySlot == 48)
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "New Piter: rho0-rhoAdv = " << m_rho0 - rhoAdv << ", pIterOld = " << pressureIterOld << ", pIterNew = " << pressureIterNew << ", deltaPIterNew-Old = " << pressureIterNew-pressureIterOld << ", vAdv=" << pTag.pointByOffset(m_vAdvOffset[m_colour]) << ", pIndex=" << i->mySlot);

	 );
      
      // --- END: new pressure value ----------------------------------------


      // --- START: iterated pair-pressure increment D_i --------------------------------------

      // recomputation with new pressure values
      computePairPressureIncrement();            

      // --- END: iterated pair-pressure increment D_i ----------------------------------------

      

      // --- START: new pressure-force increment -------------------------------------
 	 
      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(phase, m_colour, this,
	 // handles
	 const Data& pTag = i->tag;

	 // if(i->mySlot == 48)
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "PforceIncrBeforeCompute(" << i->mySlot << ")= " <<  pTag.pointByOffset(m_pforceIncrOffset));
	 
	 // dtsqF_i^Piter/m_i = 
	 // d_ii*Piter_i + sum_j [d_ij*Piter_j]
	 pTag.pointByOffset(m_pforceIncrOffset) =
	 pTag.pointByOffset(m_diiOffset)*pTag.doubleByOffset(m_pressureIterOffsetNew)
	   + pTag.pointByOffset(m_pforcePairIncrOffset);

	 // if(i->mySlot == 48)	 
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "PforceIncrAfterCompute(" << i->mySlot << ")= " <<  pTag.pointByOffset(m_pforceIncrOffset));

	 // old style (before 2017-01-18) which computed new velocity as preliminary step
	 // // v_i_adv - (d_ii*Piter_i + sum_j [d_ij*Piter_j])/dt
	 // i->v = pTag.pointByOffset(m_vAdvOffset)
	 // - (
	 //    pTag.pointByOffset(m_diiOffset)*pTag.doubleByOffset(m_pressureIterOffsetNew)
	 //    + pTag.pointByOffset(m_pforcePairIncrOffset)
	 //    )/m_dt;

	 // initialisation of new iterated density with advected density for the following summation
	 pTag.doubleByOffset(m_iterDensityOffset) = pTag.doubleByOffset(m_advDensityOffset);

	 );
	     
      // old style (before 2017-01-18) which computed new velocity as preliminary step
      // // --- END: new pressure-corrected velocity value ----------------------------------------

      // --- END: new pressure-force increment ----------------------------------------

      
      // --- START: new density due to newest pressure - self-pair-part --------
      // we must use the extra memory offset m_iterDensityOffset because we still need the
      // old density in the iteration in IntegratorIISPHconstRho::computePairPressureIncrement

      FOR_EACH_PAIR__PARALLEL
	(IntegratorIISPHconstRho, cpSelf,
	 // handles
	 const Data& p1Tag = (pair->firstPart()->tag);
	 const Data& p2Tag = (pair->secondPart()->tag);
	 const double& p1Mass = p1Tag.doubleByOffset(genMassOffsetSelf);
	 const double& p2Mass = p2Tag.doubleByOffset(genMassOffsetSelf);

	 double pairPart =
	 (p1Tag.pointByOffset(m_pforceIncrOffset) - p2Tag.pointByOffset(m_pforceIncrOffset))
	 // scalar product with the interpolation gradient
	 *-pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);

	 // ALSO DEBUG STUFF
	 double temp1 = p1Tag.doubleByOffset(m_iterDensityOffset);
	 double temp2 = p2Tag.doubleByOffset(m_iterDensityOffset);
	 
	 // if(pair->firstPart()->mySlot == 48 || pair->secondPart()->mySlot == 48)
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "selfPairContrib to density correction: p1rho-rho0Before = " << p1Tag.doubleByOffset(m_iterDensityOffset) - m_rho0 << ", p2rho-rho0Before = " << p2Tag.doubleByOffset(m_iterDensityOffset) - m_rho0);

	 // pairPart is symmetric
	 p1Tag.doubleByOffset(m_iterDensityOffset) += p2Mass*pairPart;
	 p2Tag.doubleByOffset(m_iterDensityOffset) += p1Mass*pairPart;

	 // if(pair->firstPart()->mySlot == 48 || pair->secondPart()->mySlot == 48)
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "selfPairContrib to density correction: p1("<< pair->firstPart()->mySlot << ") = " << p2Mass*pairPart <<  ", p2("<< pair->secondPart()->mySlot << ") = " << p1Mass*pairPart << ", p1rho-rho0New = " << p1Tag.doubleByOffset(m_iterDensityOffset) - m_rho0 << ", p2rho-rho0New = " << p2Tag.doubleByOffset(m_iterDensityOffset) - m_rho0 << ", rho1New-rho1Before=" << p1Tag.doubleByOffset(m_iterDensityOffset)-temp1 << ", rho2New-rho2Before=" << p2Tag.doubleByOffset(m_iterDensityOffset)-temp2);
	 
     );
      
      // --- END: new density due to newest pressure - self-pair-part --------------
      
      
      // --- START: new density due to newest pressure - non-self-pair-part --------
      
      for(vector<ColourPair*>::iterator cpit = m_mixedColourPairs.begin(); cpit != m_mixedColourPairs.end(); ++cpit)
	{
	  ColourPair* cp = (*cpit);
	  size_t genMassOffset2nd = m_genMassOffset[cp->secondColour()];
	  FOR_EACH_PAIR__PARALLEL
	    (IntegratorIISPHconstRho, cp,
	     // handles
	     const Data& p1Tag = (pair->firstPart()->tag);
	     const Data& p2Tag = (pair->secondPart()->tag);
	     const double& p2Mass = p2Tag.doubleByOffset(genMassOffset2nd);
	     
	     point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);
	     
	     p1Tag.doubleByOffset(m_iterDensityOffset) +=
	     // last product is scalar product
	     p2Mass*p1Tag.pointByOffset(m_pforceIncrOffset)*interpolGradient;
	     // no contribution to 2nd species needed so far (2017-01-11)
	     // p2Tag.doubleByOffset(m_densityOffset) +=
	     // p1Mass*p2Tag.pointByOffset(m_pforceIncrOffset[secondColour])*interpolGradient;

	     // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "nonSelfPairContrib to density correction for p1 = " << p2Mass*p1Tag.pointByOffset(m_pforceIncrOffset)*interpolGradient);
	     
	     );
	}
      
      // --- END: new density due to newest pressure - non-self-pair-part --------
      
      
      // --- START: convergence measures for density --------------------------------

      // // FOR DEBUGGING
      // size_t dbgCountBigErrs = 0;
      
      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(phase, m_colour, this,
	 double particleRelDeltaRho =
	 fabs((i->tag.doubleByOffset(m_iterDensityOffset) - m_rho0)/m_rho0);

	 // if(particleRelDeltaRho > m_avgDensityError) {
	 //   ++dbgCountBigErrs;
	 //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "particleRelDeltaRho(" << i->mySlot << ")=" << particleRelDeltaRho);
	 // }

	 deltaDensityMax = max(particleRelDeltaRho, deltaDensityMax);
	 deltaDensityAvg += particleRelDeltaRho;
	 
	 );
      
      deltaDensityAvg /= phase->returnNofPartC(m_colour);
      
      // --- END: convergence measures for density ----------------------------------

      // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "dbgCountBigErrs = " << dbgCountBigErrs);	   
      // MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "End of iteration step l = " << l << ". Here are the error limits and errors:\n"
      // 		     "maxDensityErrorLimit = " + ObjToString(m_maxDensityError) +
      // 		     " (current = " + ObjToString(deltaDensityMax) + ")\n"
      // 		     "avgDensityErrorLimit = " + ObjToString(m_avgDensityError) +
      // 		     " (current = " + ObjToString(deltaDensityAvg) + ")\n");
      
      ++l;
      
    } // end while(deltaDensityMax > m_eps...
  
  // --- END: Relaxed Jacobi iteration ----------------------------------------


  // velocity update and setting of final density and pressure
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handles
     const Data& pTag = i->tag;

     // if(i->mySlot == 48) {
     //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "velocityBeforeFinal= " <<  i->v);
     //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "vAdvBeforeFinal= " <<  pTag.pointByOffset(m_vAdvOffset[m_colour]));
     //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "PforceIncrBeforeFinal= " <<  pTag.pointByOffset(m_pforceIncrOffset));
     // }
     
     // = v_i_adv + dt*F_i^P/m_i
     // = v_i_adv - (d_ii*Piter_i + sum_j [d_ij*Piter_j])/dt
     i->v = pTag.pointByOffset(m_vAdvOffset[m_colour]) + pTag.pointByOffset(m_pforceIncrOffset)/m_dt;

     // if(i->mySlot == 48)
     //   MSG_DEBUG("IntegratorIISPHconstRho::integrateStep2", "velocityFinal= " <<  i->v);

     pTag.doubleByOffset(m_densityOffset) = pTag.doubleByOffset(m_iterDensityOffset);
     pTag.doubleByOffset(m_pressureOffset) = pTag.doubleByOffset(m_pressureIterOffsetNew);
     
     );
  
   phase->invalidate();
}


void IntegratorIISPHconstRho::integratePosition(Particle* p, Cell* cell)
{
  // size_t force_index;
  // force_index = ((Controller*) m_parent)->forceIndex();

  // const point_t& pt = p->force[force_index]/m_mass;

  // Currently (2017-01-05), wall collisions not supported by this Integrator, since
  // they would mostly violate incompressibility anyway. Therefore commented out
  // Currently (2010-05-05), pt is a const point& argument, so using it in the p->r += ... line is safe
  // cell->doCollision(p, p->r, p->v, pt, (IntegratorPosition*) this);

  // Currently (2017-01-05), standard Euler integration as in Ihmsen et al. (IEEE Transactions on Visualization and Computer Graphics 20, 426 (2014)).
  p->r += m_dt * (p->v /* + 0.5 * m_dt * pt */);

}


void IntegratorIISPHconstRho::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results)
{
  throw gError("IntegratorIISPHconstRho::solveHitTimeEquation", "Fatal ERROR: Shouldn't have been called! Contact programmer!");
}


void IntegratorIISPHconstRho::hitPos(double dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  throw gError("IntegratorIISPHconstRho::hitPos", "Fatal ERROR: Shouldn't have been called! Contact programmer!");
}


void IntegratorIISPHconstRho::computePairPressureIncrement()
{     
  ColourPair* cpSelf = M_MANAGER->cp(m_colour, m_colour);
  Phase* phase = M_PHASE;
  double dtsq = m_dt*m_dt;
  size_t genMassOffsetSelf = m_genMassOffset[m_colour];
  
  // precomputation and initialisation
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handles
     const Data& pTag = i->tag;
     const double& pDensity = pTag.doubleByOffset(m_densityOffset);

     // precomputation for next step; minus sign and dtSq already included!
     pTag.doubleByOffset(m_precomputeOffset) =
     -dtsq*pTag.doubleByOffset(genMassOffsetSelf)*pTag.doubleByOffset(m_pressureIterOffsetNew)
     /(pDensity*pDensity);

     // initialisation
     pTag.pointByOffset(m_pforcePairIncrOffset).assign(0.);

     );
  
  // No need for other ColourPairs due to currently made ASSUMPTIONS:
  // Basically we compute matrix entries for incompressibility WITHIN a SINGLE species  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorIISPHconstRho, cpSelf,
     // handles
     const Data& p1Tag = (pair->firstPart()->tag);
     const Data& p2Tag = (pair->secondPart()->tag);
     // const double& p1Mass = p1Tag.doubleByOffset(genMassOffsetSelf);
     // const double& p2Mass = p2Tag.doubleByOffset(genMassOffsetSelf);
     // const double& p1Density = p1Tag.doubleByOffset(m_densityOffset);
     // const double& p2Density = p2Tag.doubleByOffset(m_densityOffset);
     // const double& p1pressureIter = p1Tag.doubleByOffset(m_pressureIterOffsetNew);
     // const double& p2pressureIter = p2Tag.doubleByOffset(m_pressureIterOffsetNew);
     const double& p1PrecomputedNegDtSq_M_Piter_divRhoSq = p1Tag.doubleByOffset(m_precomputeOffset);
     const double& p2PrecomputedNegDtSq_M_Piter_divRhoSq = p2Tag.doubleByOffset(m_precomputeOffset);
     // the gradient of the interpolation function
     point_t interpolGradient = -pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);
     
     p1Tag.pointByOffset(m_pforcePairIncrOffset) +=
     p2PrecomputedNegDtSq_M_Piter_divRhoSq*interpolGradient;
     // ... += -dtsq*p2Mass*p2pressureIter*interpolGradient/p2Density/p2Density;
     
     // antisymmetric
     p2Tag.pointByOffset(m_pforcePairIncrOffset) -=
     p1PrecomputedNegDtSq_M_Piter_divRhoSq*interpolGradient;
     // ... -= -dtsq*p1Mass*p1pressureIter*interpolGradient/p1Density/p1Density;
     
     );
}
  
  
#ifdef _OPENMP
string IntegratorIISPHconstRho::dofIntegr() {
  return "vel_pos";
}


void IntegratorIISPHconstRho::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->force[force_index][i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
// MSG_DEBUG("IntegratorIISPHconstRho::mergeCopies", " real force to be added = " << (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] << " slot = " << p->mySlot);
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
//    MSG_DEBUG("IntegratorIISPHconstRho::mergeCopies", " force after merge = " << p->force[force_index]);
  }
}

#endif

