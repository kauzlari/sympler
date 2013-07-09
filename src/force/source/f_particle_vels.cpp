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


#include "f_particle_vels.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<FParticleVels> f_particle_vels("FParticleVels");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

//---- Constructors/Destructor ----

FParticleVels::FParticleVels(Simulation *simulation): FParticle(simulation)
{
  init();
}


FParticleVels::~FParticleVels()
{
}


void FParticleVels::computeForces(Particle* part, int force_index)
{
//   Phase *phase = ((Simulation *) m_parent)->phase();
  point_t temp;

//   FOR_EACH_FREE_PARTICLE_IN_GROUP
//     (phase, m_colour, groups,
  if (groups.empty()) {
    m_expression(&temp, part);
    part->force[force_index] += temp;

//     if (part->mySlot == 0){
//       MSG_DEBUG("FParticleVels::computeForces", "i->force = " << part->force[force_index] << "part force = " << m_is_particle_force << " pair force = " << m_is_pair_force);}

  }
  else {
    if (groups.find(part->g) != groups.end()) {
      m_expression(&temp, part);
  /*     MSG_DEBUG("FParticleVels::computeForces", "temp = " << temp);
      MSG_DEBUG("FParticleVels::computeForces", "i->tag.pointByOffset(296) = " << i->tag.pointByOffset(296));
      MSG_DEBUG("FParticleVels::computeForces", "force BEFORE = " << i->force[force_index]);*/
      part->force[force_index] += temp;
  //      MSG_DEBUG("FParticleVels::computeForces", "force AFTER = " << i->force[force_index]);
  //               );
    }
  }
}


#ifndef _OPENMP
void FParticleVels::computeForces(Pairdist* pair, int force_index)
#else
void FParticleVels::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FParticleVels::computeForces", "Fatal error: do not call FParticleVels::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FParticleVels::computeForces(int force_index)
{
  throw gError("FParticleVels::computeForces", "Fatal error: do not call FParticleVels::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FParticleVels::init()
{
  m_properties.setClassName("FParticleVels");

  m_properties.setDescription("User defined one-particle force for the momentum degree of freedom.");

  STRINGPC
    (expression, m_exprString,
     "Algebraic vector-expression for the force.");

  m_exprString = "unitVec(0)";

  m_is_pair_force = false;
  m_is_particle_force = true;
}

void FParticleVels::setup()
{
  FParticle::setup();

  m_expression.setExpression(m_exprString);
  m_expression.setReturnType(Variant::VECTOR);
  m_expression.setColour/*Pair*/(/*m_cp,*/ m_colour);
}

