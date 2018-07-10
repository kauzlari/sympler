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


#include "pair_particle_tensor.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()

const SymbolRegister<PairParticleTensor> pair_particle_tensor("PairParticleTensor");


PairParticleTensor::PairParticleTensor(Simulation* parent)
  : ValCalculatorArbitrary(parent)
{
  m_function.setReturnType(Variant::TENSOR);
  m_1stparticleFactor.setReturnType(Variant::TENSOR);
  m_2ndparticleFactor.setReturnType(Variant::TENSOR);
  m_datatype = DataFormat::TENSOR;
  init();
}

PairParticleTensor::~PairParticleTensor()
{
}

void PairParticleTensor::init()
{
//   m_properties.setClassName("PairParticleTensor");
//Not that nice but ValCalculatorArbitrary inherits from ValCalculatorPart anyways.
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("PairParticleTensor");

  m_properties.setDescription("User defined tensor property for particles, which has to be computed by pair summation. It may be defined by the attribute 'expression'. You should take care to define the attribute 'symmetry' properly. Additionally, notice that the particle expressions must be 3x3 matrices. If you only need scalars, use \"unitMat(scalar)\", where \"scalar\" is your scalar. The particle expressions are multiplied with the pair expression element by element.");

  m_1stPExpression = "unitMat(1)";
  m_2ndPExpression = "unitMat(1)";

#ifdef _OPENMP
  m_particleCalculator = true;
#endif
  m_expression = "unitMat(1)";
}

void PairParticleTensor::setup()
{
  ValCalculatorArbitrary::setup();
}

#ifdef _OPENMP
void PairParticleTensor::mergeCopies(ColourPair* cp, int thread_no) {
  size_t slot1 = m_slots.first;
  size_t slot2 = m_slots.second;

  size_t copySlot1 = m_copy_slots[thread_no].first;
  size_t copySlot2 = m_copy_slots[thread_no].second;
  size_t vecSlot1 = m_vector_slots.first;
  size_t vecSlot2 = m_vector_slots.second;

  FOR_EACH_PARTICLE_C
    (M_PHASE, cp->firstColour(),

    size_t _tmp = 0;
    for(size_t j = 0; j < SPACE_DIMS; ++j) {
      for(size_t k = 0; k < SPACE_DIMS; ++k) {
        __iSLFE->tag.tensorByOffset(slot1)(j, k) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + _tmp];
        (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + _tmp] = 0;
        ++_tmp;
      }
    }
    );

  FOR_EACH_PARTICLE_C
    (M_PHASE, cp->secondColour(),

    size_t _tmp = 0;
    for(size_t j = 0; j < SPACE_DIMS; ++j) {
      for(size_t k = 0; k < SPACE_DIMS; ++k) {
        __iSLFE->tag.tensorByOffset(slot2)(j, k) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + _tmp];
        (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + _tmp] = 0;
        ++_tmp;
      }
    }
    );
}
#endif

