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


#include "bonded_pair_particle_vector.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


const SymbolRegister<BondedPairParticleVector> bonded_pair_particle_vector("BondedPairParticleVector");


BondedPairParticleVector::BondedPairParticleVector(Simulation* parent)
  : BondedPairParticleArbitrary(parent)
{
  m_function.setReturnType(Variant::VECTOR);
  m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.setReturnType(Variant::VECTOR);
  m_datatype = DataFormat::POINT;

  init();
}

BondedPairParticleVector::~BondedPairParticleVector()
{
}

void BondedPairParticleVector::init()
{
//   m_properties.setClassName("BondedPairParticleVector");
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("BondedPairParticleVector");


  m_properties.setDescription("User defined vector property for particles, which has to be computed by summation over bonded pairs. It may be defined by the attribute 'expression'. You should take care to define the attribute 'symmetry' properly. Additionally, notice that the particle expressions must be vectors. If you only need scalars, use \"idVec(scalar)\", where \"scalar\" is your scalar. The particle expressions are multiplied with the pair expression element by element.");

  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";

// #ifdef _OPENMP
//   m_particleCalculator = true;
// #endif

  m_expression = "idVec(1)";

}

void BondedPairParticleVector::setup()
{
  BondedPairParticleArbitrary::setup(); 
}

void BondedPairParticleVector::copyMembersTo(ValCalculator* vc)
{
  BondedPairParticleArbitrary::copyMembersTo(vc);
}


#ifdef _OPENMP
void BondedPairParticleVector::mergeCopies(ColourPair* cp, int thread_no) {

  // this function currently does nothing because the parallel version computes as the serial version and writes directly in the main particle slots. The following code is from pair_particle_vector.cpp and shows how it could work

//   size_t slot1 = m_slots.first;
//   size_t slot2 = m_slots.second;

//   size_t copySlot1 = m_copy_slots[thread_no].first;
//   size_t copySlot2 = m_copy_slots[thread_no].second;
//   size_t vecSlot1 = m_vector_slots.first;
//   size_t vecSlot2 = m_vector_slots.second;

//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, cp->firstColour(),
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       i->tag.pointByOffset(slot1)[j] += (*i->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j];
//       (*i->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j] = 0;
//     }
//   );
//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, cp->secondColour(),
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       i->tag.pointByOffset(slot2)[j] += (*i->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j];
//       (*i->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j] = 0;
//     }
//   );
}
#endif
