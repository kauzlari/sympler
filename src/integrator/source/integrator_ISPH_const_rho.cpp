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


#include "integrator_ISPH_const_rho.h"

#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "controller.h"
#include "simulation.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorISPHconstRho> integrator_ISPH_const_rho("IntegratorISPHconstRho");

const point_t IntegratorISPHconstRho::dummyNullPoint = {0, 0, 0};


//---- Constructors/Destructor ----

IntegratorISPHconstRho::IntegratorISPHconstRho(Controller *controller):IntegratorPosition(controller), m_usingWalls(0), m_precomputationDone(0)
{
  init();
}


IntegratorISPHconstRho::~IntegratorISPHconstRho()
{
}


//---- Methods ----

void IntegratorISPHconstRho::init()
{
  // some modules need to know whether there is an Integrator,
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorISPHconstRho");
  
  m_properties.setDescription
    ("Incompressible SPH Integrator based on the algorithm proposed by S. Shao, E.Y.M. Lo / Advances in Water Resources 26 (2003) 787–800.\n"
     "Currently, the algorithm uses relaxed Jacobi iteration to solve the resulting pressure Poisson equation.\n"
     "NOTE THE FOLLOWING PROPERTIES, ASSUMPTIONS, LIMITATIONS, AND REQUIREMENTS:\n"
     "- This Integrator enforces incompressibility within a SINGLE species. Multi-phase incompressibility is not yet (2017-01-18) supported.\n"
     "- Walls have to be handled by introducing two additional species. With the attribute 'edgeSpecies', particles are introduced which are positioned exactly on the boundary. With the attribute 'wallSpecies', particles are introduced which are positioned \"behind\" the boundary outside of the fluid domain, usually within a shell of one maximum cutoff thickness. Particles of all species (fluid, edge, wall) must be created with ParticleCreators. The user is responsible for the correct placing of the particles. The interactions and the pressure computation between the particles are implemented as described in the reference given above, with a some differences: To achieve compatibility with arbitrarily shaped boundaries, the pressure and advected density of wall particles is computed by partial and corrected SPH-interpolation over the edge particles.\n"
     "- Hard collisions with walls or similar objects are currently (2017-12-04) not supported. Using this integrator, such collisons will not be detected and particles may penetrate walls leading to loss of particles or to their motion in 'forbidden' areas.\n"
     );

  STRINGPC
      (PPEexpression, m_PPEexpr,
       "The mathematical expression aij' for the LHS of the SPH-discretised pressure Poisson equation. The full LHS will be sumj aij(Pi-Pj) with pressure P, where aij=aij'Wij, with interpolation kernel Wij. The kernel term is multiplied automatically such that the user must only provide aij' within 'PPEexpression'.\n "
       "NOTE: Currently, aij' must be fully symmetric w.r.t. particle-index interchange i<->j."
      );
  m_PPEexpr = "undefined";
  
  STRINGPC
    (weightingFunction, m_kernelName,
     "Defines the weighting function (interpolation kernel) to be used. Note that for "
     "each pair of species, "
     "the maximum of the cutoff defined in the weighting function specified here and the cutoff used "
     "elsewhere in the simulation setup is applied.");     
  m_kernelName = "default";

  STRINGPC
    (edgeSpecies, m_edgeSpecies,
     "Name for the species of the edge particles placed directly at solid walls. Additional Integrators introducing additional degrees of freedom for the incompressible species (defined by attribute 'species') or the species defined here must be defined BELOW IntegratorISPHconstRho. Keeping the default value \"---\" will not create any species. Note that, currently, only the simultaneous definition of both the attributes 'edgeSpecies' and 'wallSpecies' is supported.");
  m_edgeSpecies = "---";

  STRINGPC
    (wallSpecies, m_wallSpecies,
     "Name for the species of the wall particles placed behind solid walls. Additional Integrators introducing additional degrees of freedom for the incompressible species (defined by attribute 'species') or the species defined here must be defined BELOW IntegratorISPHconstRho. Keeping the default value \"---\" will not create any species. Note that, currently, only the simultaneous definition of both the attributes 'edgeSpecies' and 'wallSpecies' is supported.");
  m_wallSpecies = "---";

  STRINGPC
    (edgeEdgeListName, m_edgeEdgeListName,
     "Name of the list of bonded pairs that this Integrator automatically constructs for computation of interaction terms between two edge particles."
     );
  m_edgeEdgeListName = "__ISPHedgeEdgeList";

  STRINGPC
    (edgeWallListName, m_edgeWallListName,
     "Name of the list of bonded pairs that this Integrator automatically constructs for computation of interaction terms between edge particles and wall particles."
     );
  m_edgeWallListName = "__ISPHedgeWallList";

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
    (maxErrorLimit, m_epsilonMax, 0.,
     "Convergence limit for the maximum error of the linear system "
     "solver for the PPE. The error is defined as the residual "
     "R_i = nabla.((1/rho)nabla P)_i - RHS_i for each particle i."
     );
  m_epsilonMax = 0.001;
  
  DOUBLEPC
    (avgErrorLimit, m_epsilonAvg, 0.,
     "Convergence limit for the average error of the linear system "
     "solver for the PPE. The error is defined as the residual "
     "R_i = dt^2*nabla.((1/rho)nabla P)_i - RHS_i for each particle i."
     "The average is computed by the normalised l2-norm "
     "||R||_2 = sqrt( sum_i R_i^2*dt^2 ) / N"
     );
  m_epsilonAvg = 0.0005;

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

void IntegratorISPHconstRho::setup()
{
  
  IntegratorPosition::setup();

  DataFormat::attribute_t tempAttr;

  if(m_PPEexpr == "undefined")
    throw gError("IntegratorISPHconstRho::setup", "Attribute 'PPEexpression' was not defined!");

  if(m_rho0 == HUGE_VAL)
    throw gError("IntegratorISPHconstRho::setup", "Attribute 'rho0' was not defined!");

  if(m_rho0 <= 0.)
    throw gError("IntegratorISPHconstRho::setup", "Invalid value \"" + ObjToString(m_rho0) + "\" for attribute 'rho0'. Must be >0!");

  if(m_omega <= 0.)
    throw gError("IntegratorISPHconstRho::setup", "Invalid value \"" + ObjToString(m_omega) + "\" for attribute 'omega'. Must be 0 < omega < 1!");

  if(m_omega >= 1.)
    throw gError("IntegratorISPHconstRho::setup", "Invalid value \"" + ObjToString(m_omega) + "\" for attribute 'omega'. Must be 0 < omega < 1!");

  if(m_epsilonMax <= 0.)
    throw gError("IntegratorISPHconstRho::setup", "Invalid value \"" + ObjToString(m_epsilonMax) + "\" for attribute 'maxErrorLimit'. Must be >0!");
  
  if(m_epsilonAvg <= 0.)
    throw gError("IntegratorISPHconstRho::setup", "Invalid value \"" + ObjToString(m_epsilonAvg) + "\" for attribute 'avgErrorLimit'. Must be >0!");
  
  if(m_lMax < 1)
    throw gError("IntegratorISPHconstRho::setup", "Invalid value \"" + ObjToString(m_lMax) + "\" for attribute 'nIterSteps'. Must be >0!");
  
  // interpolation kernel
  m_kernel = M_SIMULATION->findWeightingFunction(m_kernelName);
  
  // This simplifies some loops in the algorithm
  if(m_colour != 0)
    throw gError("IntegratorISPHconstRho::setup", "The incompressible species was not defined as the very first species but this is required by IntegratorISPHconstRho. To achieve this, no other modules (such as Integrators) should define a species before IntegratorISPHconstRho.");


  // -------- START: Create those other colours we interact with ----------------------

  size_t colourTester = 0;

  // YES, outside of if(m_usingWalls) since false so far and the
  // function will set m_usingWalls=true if m_edgeSpecies != "---"
  conditionalCreateColour(colourTester, m_edgeSpecies, "edgeSpecies");
  
  if(m_usingWalls)
    conditionalCreateColour(colourTester, m_wallSpecies, "wallSpecies");
  else if (m_wallSpecies != "---")
    MSG_DEBUG("IntegratorISPHconstRho::setup", "WARNING: you have set attribute 'wallSpecies' to the non-default value \"" << m_wallSpecies << "\". But this will be ignored since attribute 'edgeSpecies' was not defined.");

  if(m_usingWalls && colourTester != 2)
    throw gError("IntegratorISPHconstRho::setup", "You seem to require walls in your simulation but I was not able to create both an edge and a wall species, which are currently both required. Please check your attributes 'edgeSpecies' and 'wallSpecies' and read their documentation.");

  // -------- END: Create those other colours we interact with ----------------------
  

  // set cutoffs, activate pairs
  // at this point only this Integrator has created colours, hence fine
  // to loop over all ColourPairs
  FOR_EACH_COLOUR_PAIR
    (M_MANAGER,
     // Should not happen if m_colour == 0, right?!?
     if(cp->firstColour() != m_colour) assert(cp->secondColour() != m_colour);

     // set cutoff and activate pairs
     // might very well be that the cutoff is increased later by other modules for some ColourPairs
     cp->setCutoff(m_kernel->cutoff());
     // Here, we sitch on the non-bonded pairs for all ColourPairs.
     // Later we switch them off for edge-edge, edge-wall, wall-wall
     cp->setNeedPairs(true);
     
     );

  m_cpFluidFluid = M_MANAGER->cp(m_colour, m_colour);

  // assertion for checking preliminary limitation
  assert(m_colour == 0);

  // add internal Symbols and get offset; hence Symbols should not yet exist
  m_aiiOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, "__ISPHaii", DataFormat::DOUBLE);
  m_advDensityOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, "__ISPHadvDensity", DataFormat::DOUBLE);
  m_pressureIterOldOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, "__ISPHpIterSlot1", DataFormat::DOUBLE);
  m_pressureIterNewOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, "__ISPHpIterSlot2", DataFormat::DOUBLE);
  // We do not need the offset for m_colour. We only have it such that
  // the indexing in the m_advDensityPrecompOffset[3] array stays
  // easy
  m_advDensityPrecompOffset[m_colour] = HUGE_VAL; // not needed
  // DataFormat::addNewAttribute
  // (m_colour, "__ISPHprecompAdvDensity", DataFormat::DOUBLE);
  m_advNormalisationOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, "__ISPHnormalisation", DataFormat::DOUBLE);
  m_pressureAccelFluidOffset = DataFormat::addNewAttribute
    (m_colour, "__ISPHpressureAccel", DataFormat::POINT);
  m_pairIterContribOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, "__ISPHpairPContrib", DataFormat::DOUBLE);
  
  // pressure to be added with user defined name
  if(m_pressureName == "UNDEFINED")
    throw gError("IntegratorISPHconstRho::setup", "Attribute pressureName has value \"UNDEFINED\".");
  m_pressureOffset[m_colour] = DataFormat::addNewAttribute
    (m_colour, m_pressureName, DataFormat::DOUBLE);
  
  // the density Symbol
  if(m_densityName == "UNDEFINED")
    throw gError("IntegratorISPHconstRho::setup", "Attribute densityName has value \"UNDEFINED\".");  
  m_densityOffset = DataFormat::addNewAttribute(m_colour, m_densityName, DataFormat::DOUBLE);

  // aditional setups if using walls
  if(m_usingWalls) {
    m_cpFluidEdge = M_MANAGER->cp(m_colour, m_edgeColour);
    m_cpFluidWall = M_MANAGER->cp(m_colour, m_wallColour);
    m_cpEdgeEdge = M_MANAGER->cp(m_edgeColour, m_edgeColour);
    m_cpEdgeWall = M_MANAGER->cp(m_edgeColour, m_wallColour);
    
    // assertion for checking preliminary limitation
    assert(m_edgeColour == 1);
    assert(m_wallColour == 2);

    // add internal Symbols and get offset; hence Symbols should not yet exist
    m_aiiOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHaii", DataFormat::DOUBLE);
    m_advDensityOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHadvDensity", DataFormat::DOUBLE);
    m_advDensityOffset[m_wallColour] = DataFormat::addNewAttribute
    (m_wallColour, "__ISPHadvDensity", DataFormat::DOUBLE);
    m_pressureIterOldOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHpIterSlot1", DataFormat::DOUBLE);
    m_pressureIterNewOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHpIterSlot2", DataFormat::DOUBLE);
    m_advDensityPrecompOffset[m_edgeColour] =
      DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHprecompAdvDensity", DataFormat::DOUBLE);
    m_advDensityPrecompOffset[m_wallColour] =
      DataFormat::addNewAttribute
      (m_wallColour, "__ISPHprecompAdvDensity", DataFormat::DOUBLE);
    m_wallEdgeNormalisationOffset = DataFormat::addNewAttribute
      (m_wallColour, "__ISPHedgeNormalisation", DataFormat::DOUBLE);
    m_advNormalisationOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHnormalisation", DataFormat::DOUBLE);
    m_pairIterContribOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, "__ISPHpairPContrib", DataFormat::DOUBLE);

    // pressure to be added with user defined name
    m_pressureOffset[m_edgeColour] = DataFormat::addNewAttribute
      (m_edgeColour, m_pressureName, DataFormat::DOUBLE);
    m_pressureOffset[m_wallColour] = DataFormat::addNewAttribute
    (m_wallColour, m_pressureName, DataFormat::DOUBLE);
    
  } // end of if(m_usingWalls)
  
}


// FIXME: Can this be merged with setupAfterParticleCreation() ?
void IntegratorISPHconstRho::isAboutToStart()
{
  Phase* phase = M_PHASE;
  ManagerCell *manager = phase->manager();
  m_dt = M_CONTROLLER->dt();
  DataFormat::attribute_t tempAttr;
  
  // initialise forces to zero  
  size_t counter = 0;
  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     for (int j = 0; j < FORCE_HIST_SIZE; j++)
       __iSLFE->force[j].assign(0);
     ++counter;
    );
  if(counter == 0)
    throw gError("IntegratorISPHconstRho::isAboutToStart", "no free particles found for species " + m_species + "! Don't instantiate an Integrator for positions and velocities in that case. Use another module to create the species.");
  // FIXME: so we need some SpeciesCreator to make it more transparent
  // FIXME: put all in this function into the general setup for Nodes after the particle creation or into s.th. even more general


  if(m_usingWalls) {
    // create connected list for edge-edge
    // ColourPair* cp = 
    //   manager->cp(manager->getColour(m_edgeSpecies), manager->getColour(m_edgeSpecies));
    m_edgeEdgeListIndex = m_cpEdgeEdge->createConnectedListIndex(m_edgeEdgeListName);
    // create connected list for edge-wall
    m_edgeWallListIndex = m_cpEdgeWall->createConnectedListIndex(m_edgeWallListName);    
    
  } // end of if(m_usingWalls)
  else {
    m_edgeEdgeListIndex = HUGE_VAL;
    m_edgeWallListIndex = HUGE_VAL;
  }
  
}


void IntegratorISPHconstRho::integrateStep1()
{

  // For the following creation of bonded pairs based on non-bonded pairs, the latter ones must have been created already, which is true here, even at the first time step, when the code is in fact executed for the one and only time based on the boolean m_precomputationDone. Then, in the first call of integrateStep2(..), the true precomputations are done based on these lists

  if(m_usingWalls && !m_precomputationDone) {
  
    // loop over edge-edge in non-bonded way and create the bonded pairs
    FOR_EACH_PAIR__PARALLEL
      (IntegratorISPHconstRho, m_cpEdgeEdge,
       // handles
       Particle* p1 = pair->firstPart();
       Particle* p2 = pair->secondPart();
       m_cpEdgeEdge->addPairToConnection(p1, p2, m_edgeEdgeListIndex);     
       );

    
    // loop over edge-wall in non-bonded way and create the bonded pairs
    FOR_EACH_PAIR__PARALLEL
      (IntegratorISPHconstRho, m_cpEdgeWall,
       // handles
       Particle* p1 = pair->firstPart();
       Particle* p2 = pair->secondPart();
       m_cpEdgeWall->addPairToConnection(p1, p2, m_edgeWallListIndex);     
       );

    // Now non-bonded pairs should not be needed anymore. The following
    // line is assumed to be sufficient, such that the next call of
    // Phase::innvalidatePositions() deletes the current non-bonded
    // pairs and sets a flag which avoids that new ones are created as
    // long as needPairs remains false
    m_cpEdgeEdge->setNeedPairs(false);
    m_cpEdgeWall->setNeedPairs(false);
    // It does not hurt if the following cleaning is unnecessary.
    M_MANAGER->cp(m_wallColour, m_wallColour)->setNeedPairs(false);
    
    // Do not set m_precomputationDone = true here, because false value still required in first call of integrateStep2(..) 
    
  } // end of if(!m_precomputationDone)


  // Should not be necessary since already done after really changing
  // positions at end of integrateStep2() 
  
  // M_PHASE->invalidatePositions((IntegratorPosition*) this);

}


void IntegratorISPHconstRho::integrateStep2()
{
  PairList* pL;

  Phase *phase = M_PHASE;
  size_t force_index = M_CONTROLLER->forceIndex();
  double dtsq = m_dt*m_dt;
  double temp;
  const point_t dummyNullPoint = {0, 0, 0};

  size_t edgeAdvDensityPrecompOffset =
    m_advDensityPrecompOffset[m_edgeColour];
  size_t wallAdvDensityPrecompOffset =
    m_advDensityPrecompOffset[m_wallColour];
  size_t advDensityOffsetEdge =  m_advDensityOffset[m_edgeColour];
  size_t advDensityOffsetFluid =  m_advDensityOffset[m_colour];
  size_t advDensityOffsetWall = m_advDensityOffset[m_wallColour];
  size_t advNormalisationOffsetFluid = m_advNormalisationOffset[m_colour];
  size_t advNormalisationOffsetEdge = m_advNormalisationOffset[m_edgeColour];
  size_t aiiOffsetFluid = m_aiiOffset[m_colour];
  size_t aiiOffsetEdge = m_aiiOffset[m_edgeColour];
  
  // the following are changed during the iteration, that's why the
  // reference
  size_t& pressureIterNewFluidOffset = m_pressureIterNewOffset[m_colour];
  size_t& pressureIterOldFluidOffset = m_pressureIterOldOffset[m_colour];
  size_t& pressureIterNewEdgeOffset = m_pressureIterNewOffset[m_edgeColour];
  size_t& pressureIterOldEdgeOffset = m_pressureIterOldOffset[m_edgeColour];

  size_t pressureFluidOffset = m_pressureOffset[m_colour];
  size_t pressureEdgeOffset = m_pressureOffset[m_edgeColour];
  size_t pressureWallOffset = m_pressureOffset[m_wallColour];
  
  // works and gives the self contribution except for
  // WeightingFunctions that will complain in a clean way
  double kernelSelfContrib = m_kernel->interpolate(NULL, dummyNullPoint);
  
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handle
     const Data& pTag = i->tag;
     
     // initialisation for below
     pTag.doubleByOffset(aiiOffsetFluid) = 0.;

     // initialise advected density to self-contrib (for equal masses
     // within a species this could be much simplified, but we want to
     // keep the flexibility of different masses)
     pTag.doubleByOffset(advDensityOffsetFluid) = i->m_mass*kernelSelfContrib;
     
     // advected velocity: currently we do not need it
     // i->v += m_dt*i->force[force_index]/(pTag.doubleByOffset(genMassOffsetSelf));
     
     // advected position: commented out line is for advected v
     // i->r += m_dt*i->v;
     i->r += m_dt*(i->v + m_dt*i->force[force_index]/i->m_mass);

     );

  // Called here because positions updated here
  phase->invalidatePositions((IntegratorPosition*) this);  
  M_CONTROLLER->triggerNeighbourUpdate();

  if(m_usingWalls) {

    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_edgeColour, this,
       // handle
       const Data& pTag = i->tag;
       
       // initialisation for below
       pTag.doubleByOffset(aiiOffsetEdge) = 0.;
       
       // initialise advected density to self-contrib 
       pTag.doubleByOffset(advDensityOffsetEdge) = i->m_mass*kernelSelfContrib;
       
       // No advected velocity and position advection for m_edgeColour
       
       );
  } // end of if(m_usingWalls)
    

  // --- START: advected density based on advected positions/velocities --------------

  // fluid self-contrib already done above
  
  // advected density based on advected positions/velocities: fluid-fluid part
  advDensityPairSum(m_cpFluidFluid);

  if(m_usingWalls) {
    
    if(!m_precomputationDone) {// entered once and never again
      
      m_precomputationDone = true;  

      // density based on advected positions/velocities:
      // PRECOMPUTATIONS following
      
      // edge-edge contribs.
      advDensPrecompConnectedPairContrib(m_cpEdgeEdge, m_edgeEdgeListIndex);
      // edge-wall contribs.
      // The nominator of the contrib to wall can also be computed here
      // only the normalisation contains a dynamic variable and is
      // hence computed below, in each timestep 
      advDensPrecompConnectedPairContrib(m_cpEdgeWall, m_edgeWallListIndex);
            
      // advected density based on advected positions/velocities: edge-self-contribution 
      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(phase, m_edgeColour, this,
	 i->tag.doubleByOffset(edgeAdvDensityPrecompOffset)
	 += i->m_mass*kernelSelfContrib;
	 );
      
      // the way we compute the advected density for the wall particles does not require a self-contribution
      
    } // end of if(!m_precomputationDone)

    
    // density based on advected positions/velocities: edge precomputed contribution 
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_edgeColour, this,
       i->tag.doubleByOffset(advDensityOffsetEdge) =
       i->tag.doubleByOffset(edgeAdvDensityPrecompOffset);
       );

    // density based on advected positions/velocities: fluid-edge part;
    advDensityPairSum(m_cpFluidEdge);
          
    // advected density based on advected positions/velocities: fluid-wall part; Here, only the summation into the advected fluid density goes the usual way. The advected density of the wall particles is computed from interpolation over the edge particles;
    // contribution to fluid
    FOR_EACH_PAIR__PARALLEL
      (IntegratorISPHconstRho, m_cpFluidWall,
       
       pair->firstPart()->tag.doubleByOffset(advDensityOffsetFluid)
       += pair->secondPart()->m_mass
       *m_kernel->interpolate(pair, dummyNullPoint);

       );

    // initialise to zero for following summation
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_wallColour, this,
       i->tag.doubleByOffset(m_wallEdgeNormalisationOffset) = 0.;
       );
    
    // normalisation-denominator for wall particle variables
    // interpolated by summation over edge particles only
    // NOTE: before the following step we need the final advected
    // density for the edge particles
    pL = m_cpEdgeWall->connectedList(m_edgeWallListIndex);
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      // handles
      Particle* p1 = pair->firstPart();
      
      pair->secondPart()->tag.doubleByOffset(m_wallEdgeNormalisationOffset) +=
	p1->m_mass*m_kernel->interpolate(pair, dummyNullPoint)
	/p1->tag.doubleByOffset(advDensityOffsetEdge);
      
    }
        
    // advected density for wall particles computed from advected
    // density of edge particles (assuming Neumann=0 BC):
    // Here we only have to divide through the time-dependent
    // normalisation-denominator. The full equation is
    // rho_wall_j = sum_k(m_k*W_jk)/normalisation_j
    // where division by the previously computed
    // normalisation_j = sum_k W_jk*m_k/rho_k is done here
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_wallColour, this,
       const Data& pTag = (i->tag);
       pTag.doubleByOffset(advDensityOffsetWall) =
       pTag.doubleByOffset(wallAdvDensityPrecompOffset)
       / pTag.doubleByOffset(m_wallEdgeNormalisationOffset);
       );

  } // end of if(m_usingWalls)

  // --- END: advected density based on advected positions/velocities --------------

  
  // --- START: normalisation-denominators for fluid and edge interpolation -----

  // initialisation with fluid self-contribution
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     const Data& pTag = (i->tag);
     pTag.doubleByOffset(advNormalisationOffsetFluid) =
     i->m_mass*kernelSelfContrib
     /pTag.doubleByOffset(advDensityOffsetFluid);
     );

  // fluid-fluid
  interpolateOne(m_cpFluidFluid);

  if(m_usingWalls) {

    // initialisation with edge self-contribution
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_edgeColour, this,
       const Data& pTag = (i->tag);
       pTag.doubleByOffset(advNormalisationOffsetEdge) =
       i->m_mass*kernelSelfContrib
       /pTag.doubleByOffset(advDensityOffsetEdge);
       );
    
    // fluid-edge / edge-fluid
    interpolateOne(m_cpFluidFluid);
    
    // fluid-wall 
    FOR_EACH_PAIR__PARALLEL
      (IntegratorISPHconstRho, m_cpFluidWall,
       // handles
       Particle* p2 = pair->secondPart();
       const double& p2AdvDensity = p2->tag.doubleByOffset(advDensityOffsetEdge);
       double& p1AdvNormalisation = pair->firstPart()->tag.doubleByOffset(advNormalisationOffsetFluid);
       double kernel = m_kernel->interpolate(pair, dummyNullPoint);
       
       p1AdvNormalisation += p2->m_mass*kernel/p2AdvDensity;
       );
    
    // edge-edge
    pL = m_cpEdgeEdge->connectedList(m_edgeEdgeListIndex);
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      // handles
      Particle* p1 = pair->firstPart();
      Particle* p2 = pair->secondPart();
      const Data& p1Tag = (p1->tag);
      const Data& p2Tag = (p2->tag);
      double kernel = m_kernel->interpolate(pair, dummyNullPoint);
      
      p1Tag.doubleByOffset(advNormalisationOffsetEdge)
	+= p2->m_mass*kernel/p2Tag.doubleByOffset(advDensityOffsetEdge);
      p2Tag.doubleByOffset(advNormalisationOffsetEdge)
	+= p1->m_mass*kernel/p1Tag.doubleByOffset(advDensityOffsetEdge);
    }
    
    // edge-wall
    pL = m_cpEdgeEdge->connectedList(m_edgeWallListIndex);
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      // handles
      Particle* p2 = pair->secondPart();
      
      pair->firstPart()->tag.doubleByOffset(advNormalisationOffsetEdge)
	+= p2->m_mass*m_kernel->interpolate(pair, dummyNullPoint)
	/p2->tag.doubleByOffset(advDensityOffsetWall);
    }
    
  } // end if(m_usingWalls)
  
  // --- END: normalisation-denominators for fluid and edge interpolation -------

  
  // --- START: normalisation of advected densities for fluid and edge species --

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     const Data& pTag = (i->tag);
     pTag.doubleByOffset(advDensityOffsetFluid)
     /= pTag.doubleByOffset(advNormalisationOffsetFluid);
     );

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_edgeColour, this,
     const Data& pTag = (i->tag);
     pTag.doubleByOffset(advDensityOffsetEdge)
     /= pTag.doubleByOffset(advNormalisationOffsetEdge);
     );

  // --- END: normalisation of advected densities for fluid and edge species ----

  
  // --- START: diagonal part of a-matrix ----------------------------------------

  // For moving fluid particles: aiiFluidTot =
  // aiiFluidFluid(done) + aiiFluidEdge(done) + aiiFluidWall(done)
  // For edge particles: aiiEdgeTot = aiiEdgeFluid(done)
  // + aiiEdgeEdge(done) + aiiEdgeWall(done) + aiiEdgeWall2nd(done)
  
  // ASSUMPTION: There is only a pair-contribution and the
  // self-coontribution vanishes

  // aiiFluidFluid:
  aiiPairContrib(m_cpFluidFluid);

  if(m_usingWalls) {
    
    // aiiFluidEdge and aiiEdgeFluid
    aiiPairContrib(m_cpFluidEdge);
    
    // aiiFluidWall but NOT aiiWallFluid since wall implements boundary
    // condition for pressure and does hence not have a pressure-DOF
    FOR_EACH_PAIR__PARALLEL
      (IntegratorISPHconstRho, m_cpFluidWall,
       // handles
       const Data& p1Tag = (pair->firstPart()->tag);
       double& p1aii = p1Tag.doubleByOffset(aiiOffsetFluid);
       
       m_PPEfunc(&temp, pair);       
       p1aii += m_kernel->weight(pair, dummyNullPoint)*temp;
       );
    
    // aiiEdgeEdge
    pL = m_cpEdgeEdge->connectedList(m_edgeEdgeListIndex);
    
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {

      m_PPEfunc(&temp, pair);       
      temp *= m_kernel->weight(pair, dummyNullPoint);
      
      pair->firstPart()->tag.doubleByOffset(aiiOffsetEdge) += temp;
      pair->secondPart()->tag.doubleByOffset(aiiOffsetEdge) += temp;
    }
    
    // aiiEdgeWall: conrib to edge only
    pL = m_cpEdgeWall->connectedList(m_edgeWallListIndex);
    
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      m_PPEfunc(&temp, pair);             
      pair->firstPart()->tag.doubleByOffset(aiiOffsetEdge) +=
	m_kernel->weight(pair, dummyNullPoint) * temp;
    }
    
    // aiiEdgeWall2nd: contrib to edge only; this comes from the
    // definition of wall pressure Pj through edge pressure Pk
    // by Pj = sumk Pk*mk*Wjk/rhok/(sumk mk*Wjk/rhok).
    //      := sumk Pk*mk*Wjk/rhok/Sigmaj.
    // The resulting aii is then
    // aii = (mi/rhoi)*sumjInWall aij/Sigmaj
    //     = (mi/rhoi)*sumjInWall dWdrij*PPEexprij/Sigmaj
    pL = m_cpEdgeWall->connectedList(m_edgeWallListIndex);
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      m_PPEfunc(&temp, pair);             
      pair->firstPart()->tag.doubleByOffset(aiiOffsetEdge) +=
	// multiplication with mi/rhoi still missing here because can be
	// factored out and done right after this loop
	m_kernel->weight(pair, dummyNullPoint) * temp
	/ pair->secondPart()->tag.doubleByOffset(m_wallEdgeNormalisationOffset);
    }
    // and now the multiplication with mi/rhoi 
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_edgeColour, this,
       const Data& pTag = (i->tag);
       pTag.doubleByOffset((aiiOffsetEdge)) *=
       i->m_mass / pTag.doubleByOffset(advDensityOffsetEdge);
       );

  } // end of if(m_usingWalls)
    
  // --- END: diagonal part of a-matrix -------------------------------

  
  // --- START: The RHS -----------------------------------------------

  // The RHS used here is (rho_0-rho_adv)/rho_0/dt^2
  // So far, all the ingredients seem accessible for this RHS and no
  // further action should be necessary.
  // When trying different RHSs in different subclasses, we might at
  // least write an overloaded method that returns the RHS. 
  
  // --- END: The RHS -------------------------------------------------
  
  
  // --- START: Relaxed Jacobi iteration ------------------------------

  // initialisation of iteration
  size_t l = 0;
  bool converged = false;
  // Convergence measures. The functions computing them do the correct
  // initialisation
  double maxRes, l2Res;
       
  // initial value for pressure iteration for fluid species
  initPressure(m_colour);

  if(m_usingWalls) {

    // initial value for pressure iteration for edge species
    initPressure(m_edgeColour);

    // pressure value for wall particles by normalised interpolation
    // of pressure of edge particles
    // we compute initial value here and updated values in while loop after having new pressures
    computeWallPressure();
    
  } // end of if(m_usingWalls)  
  
  // MSG_DEBUG("IntegratorISPHconstRho::integrateStep2", "BEFORE ITER: residualEpsMax = " << residualEpsMax << ", m_maxDensityError = " << m_maxDensityError << ", residualEpsAvg = " << residualEpsAvg << ", m_avgDensityError = " << m_avgDensityError);

  // compute the aijPj i!=j contribution for initial pressure,
  // the function checks itself if only for fluid or also for edge
  totalPairContrib();

  
  // The iteration
  while(!converged || l < 2)
    {
      // MSG_DEBUG("IntegratorISPHconstRho::integrateStep2", "ITER-START, l = " << l);

      // MSG_DEBUG("IntegratorISPHconstRho::integrateStep2", "IN ITER: residualEpsMax = " << residualEpsMax << ", m_maxDensityError = " << m_maxDensityError << ", residualEpsAvg = " << residualEpsAvg << ", m_avgDensityError = " << m_avgDensityError);
      
      if(l > m_lMax) 
	throw gError("IntegratorISPHconstRho::integrateStep2", "Iteration steps exceeded user defined maximum of " + ObjToString(m_lMax) + " without convergence. Here are the error limits and errors:\n"
		     "maxResidualErrorLimit = " + ObjToString(m_epsilonMax) +
		     " (current = " + ObjToString(maxRes) + ")\n"
		     "avgResidualErrorLimit = " + ObjToString(m_epsilonAvg) +
		     " (current = " + ObjToString(l2Res) + ")\n"
		     "\nAborting.");

	
      // --- START: new pressure value ----------------------------------------

      // Set newest pressure to old; set previous old pressure to new for new summation
      size_t tempOffset = pressureIterOldFluidOffset;
      pressureIterOldFluidOffset = pressureIterNewFluidOffset;
      pressureIterNewFluidOffset = tempOffset;

      // reset to zero for new summation (we will be working on the pressureIterOffsetNew slot!!!)
      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(phase, m_colour, this,
	 // this will be used for the next pressure iteration and is hence reset here
	 i->tag.doubleByOffset(pressureIterNewFluidOffset) = 0.;
	 );

      // analogous operations in case of walls
      if(m_usingWalls) {
	
	// Set newest pressure to old; set previous old pressure to new for new summation
	tempOffset = pressureIterOldEdgeOffset;
	pressureIterOldEdgeOffset = pressureIterNewEdgeOffset;
	pressureIterNewEdgeOffset = tempOffset;
      
	// reset to zero for new summation (we will be working on the pressureIterOffsetNew slot!!!)
	FOR_EACH_FREE_PARTICLE_C__PARALLEL
	  (phase, m_edgeColour, this,
	   // this will be used for the next pressure iteration and is hence reset here
	   i->tag.doubleByOffset(pressureIterNewEdgeOffset) = 0.;
	   );

      } // end of if(m_usingWalls)

      // for fluid and edge particles: compute new pressure, i.e.
      // += RHS_i
      // *= omega/a_ii
      // += (1-omega)*Pold_i
      newPressureIter(m_colour);
      
      if(m_usingWalls) {

	newPressureIter(m_edgeColour);
	
	// recompute wall pressures based on new edge pressures
	computeWallPressure();

    } // end of if(m_usingWalls)  

      // --- END: new pressure value ----------------------------------------

      
      // Recompute aij*Pj pair contribution for residual of error
      // estimate and next iteration
      totalPairContrib();

      totalResiduals(maxRes, l2Res);

      converged = (l2Res < m_epsilonAvg && maxRes < m_epsilonMax);
      
      ++l;
      
    } // end while(residualEpsMax > m_eps...
  
  // --- END: Relaxed Jacobi iteration --------------------------------
    

  // --- START: set the final pressure --------------------------------
  
  // --- START: pressure force increment dt*FP ------------------------
  // we store this in the tag since we need it for both velocity and
  // position integration

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handles
     const Data& pTag = i->tag;

     // initialise for below
     pTag.pointByOffset(m_pressureAccelFluidOffset).assign(0.);

     // set the final pressure
     pTag.doubleByOffset(pressureFluidOffset) = pTag.doubleByOffset(pressureIterNewFluidOffset);

     );

  if(m_usingWalls) {
  
  // final pressure for edge particles
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_edgeColour, this,
     // handles
     const Data& pTag = i->tag;

     pTag.doubleByOffset(pressureEdgeOffset) = pTag.doubleByOffset(pressureIterNewEdgeOffset);

     );

  }

  // --- END: set the final pressure ----------------------------------
  
  pressureForceIncrementFluidFluid();

  if(m_usingWalls) {
    pressureForceIncrementFluidOther(m_cpFluidEdge);
    pressureForceIncrementFluidOther(m_cpFluidWall);
  }
  
  
  // --- END: pressure force increment dt*FP -------------------------


  // --- START: velocity+position update and final pressure ----------

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // handles
     const Data& pTag = i->tag;

     // velocity integration;
     // since we do not use advected velocities so far, we still have
     // to do the advection by non-pressure forces also
     i->v += m_dt*(i->force[force_index]/i->m_mass
		   + pTag.pointByOffset(m_pressureAccelFluidOffset)
		   );
  
     // The position integration we want here is
     // r(t+dt)=r(t)+(v(t)+v(t+dt)*dt/2
     // This can be rewritten as
     // r(t+dt)=rad+(FP-Fad)*dt^2/2 -------------------------> (1)  
     // with advected position rad=r(t)+dt*(v(t)+dt*Fad)
     // new pressure force FP, advection force Fad
     // and v(t+dt)=vad+dt*FP
     // with vad=v(t)+dt*Fad
     // (1) is easy to compute from what we have computed so far
     // since i->r currently holds rad
     i->r += 0.5*dtsq
     *(i->tag.pointByOffset(m_pressureAccelFluidOffset)
       - i->force[force_index]/i->m_mass // advection force
       );
     
     );
  
  // since velocities were updated
  phase->invalidateVelocities();

  // Called here because positions updated here
  phase->invalidatePositions(/* IntegratorPosition* */ this);
  M_CONTROLLER->triggerNeighbourUpdate();

  // --- END: velocity+position update and final pressure -------------
  
} // end of integrateStep2()

  
void IntegratorISPHconstRho::integratePosition(Particle* p, Cell* cell)
{
  // size_t force_index;
  // force_index = ((Controller*) m_parent)->forceIndex();

  // const point_t& pt = p->force[force_index]/m_mass;

  // Currently (2017-01-05), wall collisions not supported by this Integrator, since
  // they would mostly violate incompressibility anyway. Therefore commented out
  // Currently (2010-05-05), pt is a const point& argument, so using it in the p->r += ... line is safe
  // cell->doCollision(p, p->r, p->v, pt, (IntegratorPosition*) this);

  // The currently (2017-10-12) used algorithm does not require any position change here, since it already had to be done in the previous time step at the end of integrateStep2(). The logic of the client function Cell:updatePositions() should work fine if the call to this function makes no changes, but the particles left the calling Cell already earlier.

  // p->r += m_dt * (p->v /* + 0.5 * m_dt * pt */);

}


void IntegratorISPHconstRho::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results)
{
  throw gError("IntegratorISPHconstRho::solveHitTimeEquation", "Fatal ERROR: Shouldn't have been called! Contact programmer!");
}


void IntegratorISPHconstRho::hitPos(double dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  throw gError("IntegratorISPHconstRho::hitPos", "Fatal ERROR: Shouldn't have been called! Contact programmer!");
}

  
void IntegratorISPHconstRho::conditionalCreateColour(size_t& colourTester, string species, string attribute) {
  
  if(species == "---")
    MSG_DEBUG("IntegratorISPHconstRho::conditionalCreateColour", "User selection for "<< attribute << ": \"---\". So no species and colour will be created and considered.");
  else {
    m_usingWalls = true;
    ++colourTester;
    // adding colour and verifying that it has the right number
    assert(colourTester == Integrator::getColourAndAdd(species));
    MSG_DEBUG("IntegratorISPHconstRho::setup", "Added new species " << species << " with colour index " << colourTester << ".");
  } 
  
}


void IntegratorISPHconstRho::advDensPrecompConnectedPairContrib(ColourPair* cp, size_t listIndex) {	
  PairList* pL = cp->connectedList(listIndex);
  size_t offset1 = m_advDensityPrecompOffset[cp->firstColour()];
  size_t offset2 = m_advDensityPrecompOffset[cp->secondColour()];
  
  for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
    // handles
    Particle* p1 = pair->firstPart();
    Particle* p2 = pair->secondPart();
    double kernel = m_kernel->interpolate(pair, dummyNullPoint);
    
    p1->tag.doubleByOffset(offset1) += p2->m_mass*kernel;
    p2->tag.doubleByOffset(offset2) += p1->m_mass*kernel;
  }
  
}


void IntegratorISPHconstRho::advDensityPairSum(ColourPair* cp) {
  size_t firstAdvDensityOffset = m_advDensityOffset[cp->firstColour()];
  size_t secondAdvDensityOffset = m_advDensityOffset[cp->secondColour()];
  FOR_EACH_PAIR__PARALLEL
    (IntegratorISPHconstRho, cp,
     // handles
     Particle* p1 = pair->firstPart();
     Particle* p2 = pair->secondPart();
     double& p1AdvDensity = p1->tag.doubleByOffset(firstAdvDensityOffset);
     double& p2AdvDensity = p2->tag.doubleByOffset(secondAdvDensityOffset);
     double kernel = m_kernel->interpolate(pair, dummyNullPoint);
     
     p1AdvDensity += p2->m_mass*kernel;
     p2AdvDensity += p1->m_mass*kernel;
     );
  
}


void IntegratorISPHconstRho::interpolateOne(ColourPair* cp) {
  size_t advDensityOffset1 = m_advDensityOffset[cp->firstColour()];
  size_t advDensityOffset2 = m_advDensityOffset[cp->secondColour()];
  size_t advNormalisationOffset1 = m_advNormalisationOffset[cp->firstColour()];
  size_t advNormalisationOffset2 = m_advNormalisationOffset[cp->secondColour()];
  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorISPHconstRho, cp,
     // handles
     Particle* p1 = pair->firstPart();
     Particle* p2 = pair->secondPart();
     const Data& p1Tag = (p1->tag);
     const Data& p2Tag = (p2->tag);
     const double& p1AdvDensity = p1Tag.doubleByOffset(advDensityOffset1);
     const double& p2AdvDensity = p2Tag.doubleByOffset(advDensityOffset2);
     double& p1AdvNormalisation = p1Tag.doubleByOffset(advNormalisationOffset1);
     double& p2AdvNormalisation = p2Tag.doubleByOffset(advNormalisationOffset2);
     double kernel = m_kernel->interpolate(pair, dummyNullPoint);
     
     p1AdvNormalisation += p2->m_mass*kernel/p2AdvDensity;
     p2AdvNormalisation += p1->m_mass*kernel/p1AdvDensity;
     );
  
}


void IntegratorISPHconstRho::aiiPairContrib(ColourPair* cp) {
  
  size_t aiiOffset1 = m_aiiOffset[cp->firstColour()];
  size_t aiiOffset2 = m_aiiOffset[cp->secondColour()];
  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorISPHconstRho, cp,
     // handles
     const Data& p1Tag = (pair->firstPart()->tag);
     const Data& p2Tag = (pair->secondPart()->tag);
     double& p1aii = p1Tag.doubleByOffset(aiiOffset1);
     double& p2aii = p2Tag.doubleByOffset(aiiOffset2);
     double temp;
     m_PPEfunc(&temp, pair);
     temp *= m_kernel->weight(pair, dummyNullPoint);
     
     // ASSUMPTION: The factor within the runtime compiled expression
     // is symmetric
     p1aii += temp;
     p2aii += temp;
     );
}


void IntegratorISPHconstRho::initPressure(size_t colour) {
  Phase* phase = M_PHASE;
  size_t pressureIterOffset = m_pressureIterNewOffset[colour];
  size_t pressureOffset = m_pressureOffset[colour];
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, colour, this,
     // handles
     const Data& pTag = i->tag;
     // factor of 0.5 is best for initial value according to Ihmsen et
     // al. (IEEE Transactions on Visualization and Computer Graphics
     // 20, 426 (2014))
     pTag.doubleByOffset(pressureIterOffset)
     = 0.5*pTag.doubleByOffset(pressureOffset);
     );
}


void IntegratorISPHconstRho::computeWallPressure() {
  Phase* phase = M_PHASE;
  size_t pressureWallOffset = m_pressureOffset[m_wallColour];
  size_t pressureIterNewEdgeOffset = m_pressureIterNewOffset[m_edgeColour];
  size_t edgeAdvDensityOffset = m_advDensityOffset[m_edgeColour];
  
  // initialise to zero for following summation
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_wallColour, this,
     const Data& pTag = (i->tag);
     pTag.doubleByOffset(pressureWallOffset) = 0;
     );
  
  // non-normalised nominator
  PairList* pL = m_cpEdgeWall->connectedList(m_edgeWallListIndex);
  for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
    Particle* p1 = pair->firstPart();
    pair->secondPart()->tag.doubleByOffset(pressureWallOffset) +=
      p1->tag.doubleByOffset(pressureIterNewEdgeOffset)
      * m_kernel->interpolate(pair, dummyNullPoint) * p1->m_mass
      / p1->tag.doubleByOffset(edgeAdvDensityOffset);
  }
  
  // pressure for wall particles computed from pressure of edge
  // particles (assuming Neumann=0 BC):
  // Here we only have to divide through the time-dependent
  // normalisation-denominator. The full equation is
  // P_wall_j = sum_k(m_k*W_jk*P_k/rho_k)/normalisation_j
  // where division by the previously computed
  // normalisation_j = sum_k W_jk*m_k/rho_k is done here
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_wallColour, this,
     const Data& pTag = (i->tag);
     pTag.doubleByOffset(pressureWallOffset)
     /= pTag.doubleByOffset(m_wallEdgeNormalisationOffset);
     );
  
}


void IntegratorISPHconstRho::pairIterPContrib(ColourPair* cp) {

  size_t firstPIterNewOffset = m_pressureIterNewOffset[cp->firstColour()];
  size_t secondPIterNewOffset = m_pressureIterNewOffset[cp->secondColour()];
  size_t firstPairIterContribOffset = m_pairIterContribOffset[cp->firstColour()];
  size_t secondPairIterContribOffset = m_pairIterContribOffset[cp->secondColour()];
  // size_t firstPIterOldOffset = m_pressureIterOldOffset[cp->firstColour()];
  // size_t secondPIterOldOffset = m_pressureIterOldOffset[cp->secondColour()];
  double temp;

  FOR_EACH_PAIR__PARALLEL
    (IntegratorISPHconstRho, cp,
     // handles
     const Data& p1Tag = (pair->firstPart()->tag);
     const Data& p2Tag = (pair->secondPart()->tag);
     const double& p1PIterNew = p1Tag.doubleByOffset(firstPIterNewOffset);
     const double& p2PIterNew = p2Tag.doubleByOffset(secondPIterNewOffset);
     double& p1PairIterContrib = p1Tag.doubleByOffset(firstPairIterContribOffset);
     double& p2PairIterContrib = p2Tag.doubleByOffset(secondPairIterContribOffset);
     // const double& p1PIterOld = p1Tag.doubleByOffset(firstPIterOldOffset);
     // const double& p2PIterOld = p2Tag.doubleByOffset(secondPIterOldOffset);
     m_PPEfunc(&temp, pair);
     temp *= m_kernel->weight(pair, dummyNullPoint);

     p1PairIterContrib -=  temp*p2PIterNew;
     p2PairIterContrib -=  temp*p1PIterNew;
     );
  
}


/*!
 * Computes the final new pressure for the given iteration step:
 * Adds the RHS (\a m_rho0 - rho_adv) / (\a m_rho0 * dt^2),  
 * multiplies with \a m_omega / aii
 * finally adds (1-omega)*Pold.
 * ASSUMPTION: The operations are performed on the data in the 
 * \a Particle tag accessed by \a m_pressureIterOldOffset for \a colour.
 * @param colour Colour of the particles the function should operate on
 */
void IntegratorISPHconstRho::newPressureIter(size_t colour) {

  Phase* phase = M_PHASE;
  double temp = m_omega/m_dt/m_dt/m_rho0;
  double oneMinOm = 1-m_omega;
  size_t pressureIterOldOffset = m_pressureIterOldOffset[colour];
  size_t pressureIterNewOffset = m_pressureIterNewOffset[colour];
  size_t pairContribOffset = m_pairIterContribOffset[colour];
  size_t advDensityOffset = m_advDensityOffset[colour];
  size_t aiiOffset = m_aiiOffset[colour];

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, colour, this,
     const Data& pTag = i->tag;
     double& newP = pTag.doubleByOffset(pressureIterNewOffset);
     // add pair contrib + RHS without denominator
     newP += pTag.doubleByOffset(pairContribOffset)
     + m_rho0 - pTag.doubleByOffset(advDensityOffset);
     // *=omega/(dt^2*rho0)
     newP *= temp/pTag.doubleByOffset(aiiOffset);
     // += (1-omega)*Pold
     newP += oneMinOm*pTag.doubleByOffset(pressureIterOldOffset);
     );
}


void IntegratorISPHconstRho::totalPairContrib() {
  // ASSUMPTION: No self-contribution terms needed because each
  // sum_j=i(aij=i*Pj=i) is zero because of aii=0
  // And, yes, there is a aii*Pi term, where aii=sumj aij !=0, but
  // that's a different story and already done above

  Phase* phase = M_PHASE;
  size_t pairContribFluidOffset = m_pairIterContribOffset[m_colour];
  size_t pairContribEdgeOffset = m_pairIterContribOffset[m_edgeColour];
  size_t pressureIterNewEdgeOffset = m_pressureIterNewOffset[m_edgeColour];
  size_t advDensityEdgeOffset = m_advDensityOffset[m_edgeColour];
  size_t pressureWallOffset = m_pressureOffset[m_wallColour];
  double temp;
  
  // reset to zero for new summation (we will be working on the pressureIterOffsetNew slot!!!)
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // this will be used for the next pressure iteration and is hence reset here
     i->tag.doubleByOffset(pairContribFluidOffset) = 0.;
     );
  
  // fluid-fluid pressure contribution
  pairIterPContrib(m_cpFluidFluid);
  
  // analogous operations in case of walls
  if(m_usingWalls) {
	
    // reset to zero for new summation (we will be working on the pressureIterOffsetNew slot!!!)
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_edgeColour, this,
       // this will be used for the next pressure iteration and is hence reset here
       i->tag.doubleByOffset(pairContribEdgeOffset) = 0.;
       );
    
    // subtraction term for edge-wall pressure contribution
    // sum_j a_ij*(1/Sigma_j)*(-1)*(m_i/rho_i)*W_ji*P_i^l
    // with i edge and j wall
    // This is done first for edge since we want to multiply the
    // i-terms in a particle loop which can only be done, if this
    // is the only contribution so far
    PairList* pL = m_cpEdgeWall->connectedList(m_edgeWallListIndex);	
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      m_PPEfunc(&temp, pair);
      pair->firstPart()->tag.doubleByOffset(pairContribEdgeOffset)
	// yes, +=, since we substract from terms that are all
	// added by -=
	// this is Wji = Wij
	+= m_kernel->interpolate(pair, dummyNullPoint)
	// this is aij
	* temp * m_kernel->weight(pair, dummyNullPoint);
    }
    // now multiply the particle-prefactor (m_i/rho_i)*P_i^l
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_edgeColour, this,
       const Data& pTag = i->tag;
       pTag.doubleByOffset(pairContribEdgeOffset)
       *= i->m_mass
       * pTag.doubleByOffset(pressureIterNewEdgeOffset)
       / pTag.doubleByOffset(advDensityEdgeOffset);
       );

    // fluid-edge / edge-fluid pressure contribution
    pairIterPContrib(m_cpFluidEdge);      
      
    // fluid-wall pressure contribution
    FOR_EACH_PAIR__PARALLEL
      (IntegratorISPHconstRho, m_cpFluidWall,

       m_PPEfunc(&temp, pair);

       pair->firstPart()->tag.doubleByOffset(pairContribFluidOffset)
       -=  m_kernel->weight(pair, dummyNullPoint) * temp
       * pair->secondPart()->tag.doubleByOffset(pressureWallOffset);
       
       );

    // edge-edge pressure contribution
    pL = m_cpEdgeEdge->connectedList(m_edgeEdgeListIndex);
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      // handles
      const Data& p1Tag = (pair->firstPart()->tag);
      const Data& p2Tag = (pair->secondPart()->tag);
      
      m_PPEfunc(&temp, pair);
      temp *= m_kernel->weight(pair, dummyNullPoint);
	  
      p1Tag.doubleByOffset(pairContribEdgeOffset) -= temp*p2Tag.doubleByOffset(pressureIterNewEdgeOffset);
      p2Tag.doubleByOffset(pairContribEdgeOffset) -= temp*p1Tag.doubleByOffset(pressureIterNewEdgeOffset);
    }
	
    // edge-wall pressure contribution
    // sum_j a_ij*(1/Sigma_j)*sum_k (m_k/rho_k)*W_jk*P_k^l
    // (with j in wall and i,k in edge)
    // = sum_j a_ij*P_j^l
    // where P_j^l was precomputed for this l
    // NOTE: The subtraction necessary for this term was already
    // done above.
    pL = m_cpEdgeWall->connectedList(m_edgeWallListIndex);	
    for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
      m_PPEfunc(&temp, pair);
      pair->firstPart()->tag.doubleByOffset(pairContribEdgeOffset)
	-= m_kernel->weight(pair, dummyNullPoint) * temp
	* pair->secondPart()->tag.doubleByOffset(pressureWallOffset);
    }
    
  } // end of if(m_usingWalls)
  
} // end of void IntegratorISPHconstRho::totalPairContrib()

/*!
 * Computes total normalised l2-norm and the infinity norm    
 */
void IntegratorISPHconstRho::totalResiduals(double& maxRes, double& l2Res) {

  Phase* phase = M_PHASE;

  size_t contribs = phase->returnNofPartC(m_colour);
  
  maxRes = 0.;
  l2Res = 0.;

  residualsPerColour(m_colour, maxRes, l2Res);
    
  if(m_usingWalls) {

    contribs += phase->returnNofPartC(m_edgeColour);

    residualsPerColour(m_edgeColour, maxRes, l2Res);
    
  } // end of if(m_usingWalls)

  // normalised l2-norm
  l2Res = sqrt(l2Res)/contribs;
  
}


/*!
 * Compute maximum and squared l2-norm for residuals R_i`=Ri*dt^2 with  
 * R_i=RHS_i-AijPj = RHS_i - (AijPjpairs + AiiPi) for one colour  
 */
void IntegratorISPHconstRho::residualsPerColour(size_t colour, double& maxRes, double& sqL2Res) {

  Phase* phase = M_PHASE;
  double temp;
  double dtSq = m_dt*m_dt;
  double rhsDenominator = dtSq*m_rho0; 
  // size_t rhsOffset = m_rhsOffset[colour];
  size_t pairIterContribOffset = m_pairIterContribOffset[colour];
  size_t advDensityOffset = m_advDensityOffset[colour];
  size_t aiiOffset = m_aiiOffset[colour];
  size_t pressureIterNewOffset = m_pressureIterNewOffset[colour];
  
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, colour, this,
     const Data& pTag = i->tag;
     
     temp =
     // RHS
     // pTag.doubleByOffset(rhsOffset)
     (m_rho0 - pTag.doubleByOffset(advDensityOffset))/rhsDenominator
     - (// Aij*Pj
	pTag.doubleByOffset(pairIterContribOffset)
	// Aii*Pi
	+ pTag.doubleByOffset(aiiOffset)
	*pTag.doubleByOffset(pressureIterNewOffset)
	);

     temp*=dtSq;
     
     maxRes = max(fabs(temp), maxRes);

     temp*=temp;     
     sqL2Res += temp;
     
     );
  
}


/*!
 * Pressure force acceleration contribution FP/m of particles from 
 * other colours to the fluid particles.
 * We use the SPH-discretisation 
 * dv/dt = -nablaP/rho = -sum_j[mj(Pi/rhoi^2+Pj/rhoj^2)]nablaWij
 * The density inserted for rhoi is currently (2017-12-01) the advected 
 * density.
 */
void IntegratorISPHconstRho::pressureForceIncrementFluidOther(ColourPair* cp) {

  size_t advDensityFluidOffset = m_advDensityOffset[m_colour];
  size_t advDensityOtherOffset = m_advDensityOffset[cp->secondColour()];
  size_t pressureFluidOffset = m_pressureOffset[m_colour];
  size_t pressureOtherOffset = m_pressureOffset[cp->secondColour()];
  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorISPHconstRho, cp,
     Particle* p1 = pair->firstPart();
     Particle* p2 = pair->secondPart();
     const Data& p1Tag = p1->tag;
     const Data& p2Tag = p2->tag;
     const double& p1Density = p1Tag.doubleByOffset(advDensityFluidOffset);
     const double& p2Density = p2Tag.doubleByOffset(advDensityOtherOffset);
     
     p1Tag.pointByOffset(m_pressureAccelFluidOffset)
     += p2->m_mass*
     (p1Tag.doubleByOffset(pressureFluidOffset)/p1Density/p1Density
      + p2Tag.doubleByOffset(pressureOtherOffset)/p2Density/p2Density)
     // negative kernel gradient
     * pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);
     
     );
  
}


/*!
 * Pressure force acceleration contribution dt*FP/m of fluid particles 
 * to the fluid particles.
 * We use the SPH-discretisation 
 * dv/dt = -nablaP/rho = -sum_j[mj(Pi/rhoi^2+Pj/rhoj^2)]nablaWij
 * The density inserted for rhoi is currently (2017-12-01) the advected 
 * density.
 */
void IntegratorISPHconstRho::pressureForceIncrementFluidFluid() {

  size_t advDensityFluidOffset = m_advDensityOffset[m_colour];
  size_t pressureFluidOffset = m_pressureOffset[m_colour];
  
  FOR_EACH_PAIR__PARALLEL
    (IntegratorISPHconstRho, m_cpFluidFluid,
     Particle* p1 = pair->firstPart();
     Particle* p2 = pair->secondPart();
     const Data& p1Tag = p1->tag;
     const Data& p2Tag = p2->tag;
     const double& p1Density = p1Tag.doubleByOffset(advDensityFluidOffset);
     const double& p2Density = p2Tag.doubleByOffset(advDensityFluidOffset);
     
     // note that the m_dt is already included
     point_t temp =
     (p1Tag.doubleByOffset(pressureFluidOffset)/p1Density/p1Density
      + p2Tag.doubleByOffset(pressureFluidOffset)/p2Density/p2Density)
     // negative kernel gradient
     * pair->cartesian()*m_kernel->weight(pair, dummyNullPoint);
     
     p1Tag.pointByOffset(m_pressureAccelFluidOffset) += p2->m_mass*temp;
     // antisymmetric
     p2Tag.pointByOffset(m_pressureAccelFluidOffset) -= p1->m_mass*temp;
     
     ); 
  
}


#ifdef _OPENMP
string IntegratorISPHconstRho::dofIntegr() {
  return "vel_pos";
}


void IntegratorISPHconstRho::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->force[force_index][i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
// MSG_DEBUG("IntegratorISPHconstRho::mergeCopies", " real force to be added = " << (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] << " slot = " << p->mySlot);
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
//    MSG_DEBUG("IntegratorISPHconstRho::mergeCopies", " force after merge = " << p->force[force_index]);
  }
}

#endif

