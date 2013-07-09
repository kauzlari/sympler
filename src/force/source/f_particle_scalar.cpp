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


#include "f_particle_scalar.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<FParticleScalar> f_particle_scalar("FParticleScalar");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

//---- Constructors/Destructor ----

FParticleScalar::FParticleScalar(Simulation *simulation): FParticle(simulation)
{
  init();
} 


FParticleScalar::~FParticleScalar()
{
}


void FParticleScalar::computeForces(Particle* part, int force_index)
{
  Phase *phase = ((Simulation *) m_parent)->phase();

  // necessary, because force could depend on pair-properties
//   M_PAIRCREATOR->createDistances();

  double temp;

//   FOR_EACH_FREE_PARTICLE_IN_GROUP
//       (phase, m_colour, groups,

  if (groups.empty()) {
    m_expression(&temp, part);
//     MSG_DEBUG("FParticleScalar::computeForces", "force BEFORE = " << part->tag.doubleByOffset(m_force_offset[force_index]));
    part->tag.doubleByOffset(m_force_offset[force_index]) += temp;

//     MSG_DEBUG("FParticleScalar::computeForces", "force AFTER = " << part->tag.doubleByOffset(m_force_offset[force_index]));
  }
  else {
    if (groups.find(part->g) != groups.end()) {
        m_expression(&temp, part);
  /*     MSG_DEBUG("FParticleScalar::computeForces", "temp = " << temp);
        MSG_DEBUG("FParticleScalar::computeForces", "i->tag.pointByOffset(296) = " << i->tag.pointByOffset(296));
        MSG_DEBUG("FParticleScalar::computeForces", "force BEFORE = " << i->force[force_index]);*/
        part->tag.doubleByOffset(m_force_offset[force_index]) += temp;
//         MSG_DEBUG("FParticleScalar::computeForces", "force AFTER = " << part->tag.doubleByOffset(m_force_offset[force_index]));
  //       );
    }
  }
}


#ifndef _OPENMP
void FParticleScalar::computeForces(Pairdist* pair, int force_index)
#else
void FParticleScalar::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FParticleScalar::computeForces", "Fatal error: do not call FParticleScalar::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FParticleScalar::computeForces(int force_index)
{
  throw gError("FParticleScalar::computeForces", "Fatal error: do not call FParticleScalar::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FParticleScalar::init()
{
  m_properties.setClassName("FParticleScalar");

  m_properties.setDescription("User defined one-particle force for a scalar degree of freedom.");

  STRINGPC
      (expression, m_exprString,
       "Algebraic scalar expression for the force.");
  STRINGPC
      (scalar, m_scalar_name,
       "Name of the scalar field this force acts on.");

  m_scalar_name = "undefined";

  m_exprString = "0";
  m_is_pair_force = false;
  m_is_particle_force = true;
}

void FParticleScalar::setup()
{
  FParticle::setup();

  if(m_scalar_name == "undefined")
    throw gError("FParticleScalar::setup", "Attribute 'scalar' has value \"undefined\".");


  m_expression.setExpression(m_exprString);
  m_expression.setReturnType(Variant::SCALAR);
  m_expression.setColour/*Pair*/(/*m_cp,*/ m_colour);

  DataFormat::attribute_t attr =
      Particle::s_tag_format[m_colour].attrByName(m_scalar_name);


  if(attr.datatype != DataFormat::DOUBLE)
    throw gError("FParticleScalar::setup", "the symbol " + m_scalar_name +
        " is registerd as a non-scalar for species " + m_species + "."
                );

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i] =
        Particle::s_tag_format[m_colour].attrByName(string("force_"
        + m_scalar_name + "_" + ObjToString(i))).offset;
  }

}

