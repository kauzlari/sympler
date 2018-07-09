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


#include "f_particle_tensor.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<FParticleTensor> f_particle_tensor("FParticleTensor");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

//---- Constructors/Destructor ----

FParticleTensor::FParticleTensor(Simulation *simulation): FParticle(simulation)
{
  init();
} 


FParticleTensor::~FParticleTensor()
{
}


void FParticleTensor::computeForces(Particle* part, int force_index)
{
  Phase *phase = ((Simulation *) m_parent)->phase();

  // necessary, because force could depend on pair-properties 
//   M_PAIRCREATOR->createDistances();
  
  tensor_t temp;

//   FOR_EACH_FREE_PARTICLE_IN_GROUP
//       (phase, m_colour, groups,

  if (groups.empty()) {
    m_expression(&temp, part);
    part->tag.tensorByOffset(m_force_offset[force_index]) += temp;
//MSG_DEBUG("FParticleTensor::computeForces", "force AFTER = " << part->tag.tensorByOffset(m_force_offset[force_index]));
  }
  else {
    if (groups.find(part->g) != groups.end()) {
        m_expression(&temp, part);
  /*     MSG_DEBUG("FParticleTensor::computeForces", "temp = " << temp);
        MSG_DEBUG("FParticleTensor::computeForces", "i->tag.pointByOffset(296) = " << i->tag.pointByOffset(296));
        MSG_DEBUG("FParticleTensor::computeForces", "force BEFORE = " << i->force[force_index]);*/
        part->tag.tensorByOffset(m_force_offset[force_index]) += temp;
        
      //  MSG_DEBUG("FParticleTensor::computeForces", "force AFTER = " << part->tag.tensorByOffset(m_force_offset[force_index]));
  //       );
    }
  }
}


#ifndef _OPENMP
void FParticleTensor::computeForces(Pairdist* pair, int force_index)
#else
void FParticleTensor::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FParticleTensor::computeForces", "Fatal error: do not call FParticleTensor::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");  
}


void FParticleTensor::computeForces(int force_index)
{
  throw gError("FParticleTensor::computeForces", "Fatal error: do not call FParticleTensor::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");  
}


void FParticleTensor::init()
{
  m_properties.setClassName("FParticleTensor");

  m_properties.setDescription("User defined one-particle force for a tensor degree of freedom.");

  STRINGPC
      (expression, m_exprString,
       "Algebraic tensor expression for the force.");
  STRINGPC
      (tensor, m_tensor_name,
       "Name of the tensor field this force acts on.");

  m_tensor_name = "undefined";
  
  m_exprString = "unitMat(0)";
  m_is_pair_force = false;
  m_is_particle_force = true;
}

void FParticleTensor::setup()
{
  FParticle::setup();

  if(m_tensor_name == "undefined")
    throw gError("FParticleTensor::setup", "Attribute 'tensor' has value \"undefined\".");
    
  
  m_expression.setExpression(m_exprString);
  m_expression.setReturnType(Variant::TENSOR);
  m_expression.setColour/*Pair*/(/*m_cp,*/ m_colour);
  
  DataFormat::attribute_t attr =
      Particle::s_tag_format[m_colour].attrByName(m_tensor_name);
  
  
  if(attr.datatype != DataFormat::TENSOR) 
    throw gError("FParticleTensor::setup", "the symbol " + m_tensor_name + 
        " is registerd as a non-tensor for species " + m_species + "."
                );

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i] =
        Particle::s_tag_format[m_colour].attrByName(string("force_" 
        + m_tensor_name + "_" + ObjToString(i))).offset;
  }      
 
}

