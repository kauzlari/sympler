/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
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



#include "lennard_jones.h"

#include "phase.h"
#include "simulation.h"
#include "pair_creator.h"
#include "threads.h"


using namespace std;

const GenFTypeConcr<LJ> lj("LJ");


#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()


//---- Constructors/Destructor ----

LJ::LJ(Simulation *simulation)
    : FPair(simulation)
{
    init();
}


// void LJ::computeEnergy(double& energy, group_t *for_groups)
// {   
// //   Phase *phase = M_SIMULATION->phase();
//   double temp = 0;
//  M_PAIRCREATOR->createDistances();

//   ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));    
//   FOR_EACH_PAIR_HALF_IN_GROUP
//     (m_cp, (*for_groups),
//      double r = pair->abs();
//      if (r < m_cutoff) {
//        double r6i = 1/(r*r);
//        r6i = r6i*r6i*r6i;
// //        double r6i = pair->tag.doubleByOffset(m_compute_r6i_offset);
//        temp += pair_factor * (4*m_epsilon*(m_sigma_pow_12*r6i*r6i-m_sigma_pow_6*r6i) - m_shift_energy);
//      }
//     );
//   temp *= 4;
//   energy += temp;
// }


void LJ::init()
{
  m_properties.setClassName("LJ");

  m_properties.setDescription(
    "Force based on Lennard-Jones potential. Parametrization is "
    "chosen to be: V(r)=4*epsilon*[(sigma/r)^12-(sigma/r)^6]"
  );

  DOUBLEPC
    (sigma, m_sigma, 0,
     "Sigma parameter of the Lennard-Jones potential. Should be >0.");

  DOUBLEPC
    (epsilon, m_epsilon, 0,
     "Epsilon parameter of the Lennard-Jones potential. Should be >0.");

//   BOOLPC
//     (shiftPotential, m_shift_potential,
//      "Shift the potential to drop to zero smoothly at the cutoff radius.");

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Cutoff distance. Should be >0.");

  m_sigma = 1;
  m_epsilon = 1;
  m_cutoff = 2.5;
//   m_shift_potential = true;
}


void LJ::setup()
{
  FPair::setup();
  ColourPair *cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_half_sigma_pow_6 = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma;
  m_sigma_pow_12 = m_half_sigma_pow_6*m_half_sigma_pow_6;
  // that's how it is needed in the force
  m_half_sigma_pow_6 *= 0.5;
  
//   if (m_shift_potential) {
//     double rc6i = 1 / (m_cutoff*m_cutoff*m_cutoff*m_cutoff*m_cutoff*m_cutoff);
        
//     m_shift_energy = 4*m_epsilon*(m_sigma_pow_12*rc6i*rc6i-m_sigma_pow_6*rc6i);
//   }
  
  cp->setCutoff(m_cutoff);

  m_c1 = M_MANAGER->getColour(m_species.first);
  m_c2 = M_MANAGER->getColour(m_species.second);

}

void LJ::computeForces(Particle* part, int force_index)
{
  throw gError("LJ::computeForces", "Fatal error: do not call FPairVels::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void LJ::computeForces(int force_index)
{
  throw gError("LJ::computeForces", "Fatal error: do not call FPairVels::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}



#ifdef _OPENMP
void LJ::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof = "vel_pos";
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
      if (col1 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();
      }
      else if (col1 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();
      }
      else {
        throw gError("LJ::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
  if (col2 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
      if (col2 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();
      }
      else if (col2 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();        
      }
      else {
        throw gError("LJ::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
}

#endif
