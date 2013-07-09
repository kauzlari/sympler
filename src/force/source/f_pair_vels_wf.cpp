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



#include "f_pair_vels_wf.h"

#include "threads.h"
#include "particle.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_PAIRCREATOR M_PHASE->pairCreator()
#define M_MANAGER M_PHASE->manager()

const GenFTypeConcr<FPairVelsWF> f_pair_vels_wf("FPairVelsWF");

//---- Constructors/Destructor ----

FPairVelsWF::FPairVelsWF(Simulation *simulation): FPairArbitraryWF(simulation)
{
  init();
}


FPairVelsWF::~FPairVelsWF()
{
}


void FPairVelsWF::init()
{
  m_properties.setClassName("FPairVelsWF");

  m_properties.setDescription(
    "This is a completely general pair force FPV on a velocity-vector v_i such that:\n"
 "  dv_i = FPV*dt\n"
 "       = Sum_j(particleFactor(i)_ij*pairFactor_ij*weight_ij)*dt\n"
 "where pairFactor_ij includes all (anti-) symmetric pair contributions of the pair ij,"
 " and particleFactor(i)_ij the non-symmetric ones for particle i."
 " weight_ij represents the derivative of the used weighting function."
 " Note that the expression has to give a vector in the end,"
 " otherwise the compiler will report an error.");

  m_symmetry = -1;

  m_pairFactorStr = "idVec(1)";
  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";

  m_is_particle_force = false;
  m_is_pair_force = true;
}


//---- Methods ----

void FPairVelsWF::computeForces(Particle* part, int force_index)
{
  throw gError("FPairVelsWF::computeForces", "Fatal error: do not call FPairVelsWF::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairVelsWF::computeForces(int force_index)
{
  throw gError("FPairVelsWF::computeForces", "Fatal error: do not call FPairVelsWF::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairVelsWF::setup()
{
  FPairArbitraryWF::setup();

  m_pairFactor.setReturnType(Variant::VECTOR);
  m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.setReturnType(Variant::VECTOR);

// Setting the colours this force works on
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if ((col1 == cp->secondColour()) && (col2 == cp->firstColour())) {
    size_t dummy = col1;
    col1 = col2;
    col2 = dummy;
  }
  else if ((col1 == cp->firstColour()) && (col2 == cp->secondColour())) {
//    MSG_DEBUG("FPairVelsWF::setup", "Force colours same order as ColourPair's colours");
  }
  else {
    throw gError("FPairVelsWF::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FPairVelsWF::setForceSlots(Integrator* intr, int thread_no) {
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
        throw gError("FPairVelsWF::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
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
        throw gError("FPairVelsWF::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
}

#endif
