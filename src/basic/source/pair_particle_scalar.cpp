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


#include "pair_particle_scalar.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()

const SymbolRegister<PairParticleScalar> pair_particle_scalar("PairParticleScalar");


PairParticleScalar::PairParticleScalar(Simulation* parent)
  : ValCalculatorArbitrary(parent)
{
  m_function.setReturnType(Variant::SCALAR);
  m_1stparticleFactor.setReturnType(Variant::SCALAR);
  m_2ndparticleFactor.setReturnType(Variant::SCALAR);
  m_datatype = DataFormat::DOUBLE;
  init();
}

PairParticleScalar::~PairParticleScalar()
{
}

void PairParticleScalar::init()
{
//   m_properties.setClassName("PairParticleScalar");
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("PairParticleScalar");

  m_properties.setDescription("User defined scalar property for particles, which has to be computed by pair summation. It may be defined by the attribute 'expression'. Additionally you should take care to define the attribute 'symmetry' properly.");

  m_1stPExpression = "1";
  m_2ndPExpression = "1";

#ifdef _OPENMP
  m_particleCalculator = true;
#endif
  m_expression = "1";
}

void PairParticleScalar::setup()
{
  ValCalculatorArbitrary::setup();
}

#ifdef _OPENMP
void PairParticleScalar::mergeCopies(ColourPair* cp, int thread_no) {
  size_t slot1 = m_slots.first;
  size_t slot2 = m_slots.second;

  size_t copySlot1 = m_copy_slots[thread_no].first;
  size_t copySlot2 = m_copy_slots[thread_no].second;
  size_t vecSlot1 = m_vector_slots.first;
  size_t vecSlot2 = m_vector_slots.second;

  FOR_EACH_PARTICLE_C
  (M_PHASE, cp->firstColour(),
    __iSLFE->tag.doubleByOffset(slot1) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1];
      (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1] = 0;
  );
  FOR_EACH_PARTICLE_C
  (M_PHASE, cp->secondColour(),
    __iSLFE->tag.doubleByOffset(slot2) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2];
      (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2] = 0;
  );
}
#endif

