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


#include "bonded_pair_particle_tensor_noise_vector.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


const SymbolRegister<BondedPairParticleTensorNoiseVector> bonded_pair_particle_tensor_noise_vector("BondedPairParticleTensorNoiseVector");


BondedPairParticleTensorNoiseVector::BondedPairParticleTensorNoiseVector(Simulation* parent)
  : BondedPairParticleCalc(parent), m_controller(NULL)
{

  m_datatype = DataFormat::POINT;
  m_stage = 0;
  init();
}

BondedPairParticleTensorNoiseVector::~BondedPairParticleTensorNoiseVector()
{
}

void BondedPairParticleTensorNoiseVector::init()
{
  //   m_properties.setClassName("BondedPairParticleTensorNoiseVector");
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("BondedPairParticleTensorNoiseVector");
  
  
  m_properties.setDescription
    (
     "Calculator doing summation over bonded pairs and computing a function of a matrix \"{dW}\" of independent Wiener increments. The summation increment for particle i is \n"     
     "dx_i = particleFactor_i({dW})*Sum_j(pairFactor_ij({dW}))*dt\n"
     "where pairFactor_ij({dW}) includes all pair contributions of the pair ij,\n"
     "particleFactor_i(j)({dW}) represents factors specific to particle i(j),\n"
     "and {dW} is the symbol for the matrix Wiener increment, which may be used "
     "in the expressions as indicated. Further note that this module, as a preparation for Velocity-Verlet time integration, divides its contribution by sqrt(dt), with dt the integration-timestep as defined in the Controller."
);

  STRINGPC  
      (pairFactor, m_pairFactorStr,
       "Function for pairFactor_ij. Type some nonsense to \n"
       "obtain a complete list of possible variables and constants.\n"
       "The expression may contain scalars, vectors and tensors,"
       " but as a whole it must represent the type of the corresponding force.");
  
  STRINGPC
      (particleFactor_i, m_1stPExpression,
       "The mathematical expression of the additional particle factor for the first particle."
      );

  STRINGPC
      (particleFactor_j, m_2ndPExpression,
       "The mathematical expression of the additional particle factor for the second particle."
      );

  m_pairFactorStr = "idVec(1)";
  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";

  
//   DOUBLEPC
//     (noise, m_noise, 0,
//      "Noise amplitude.");
  
//   m_noise = 1;
  
  BOOLPC(randomize, m_randomize, "Initalize random number generator with a non-default seed. This makes the outcome of two simulations with identical setup different.");
  
  m_randomize = false;

}

void BondedPairParticleTensorNoiseVector::setup()
{
  MSG_DEBUG("BondedPairParticleTensorNoiseVector::setup", ": START: m_symmetry = " << m_symmetry << ", &m_symmetry=" << &m_symmetry);    

  m_controller = M_CONTROLLER;

//   m_noise_and_time = m_noise/sqrt(M_SIMULATION->controller()->dt());
  

  ColourPair *cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  assert(cp);
  
  m_pairFactor.setExpression(m_pairFactorStr);
  m_pairFactor.setColourPair(cp);
  m_1stparticleFactor.setExpression(m_1stPExpression);
  m_1stparticleFactor.setColourPair(cp);
  m_2ndparticleFactor.setExpression(m_2ndPExpression);
  m_2ndparticleFactor.setColourPair(cp);
  if (m_1stparticleFactor.isNull()) {
    throw gError("BondedPairParticleTensorNoiseVector::setup", "For " + m_properties.className() + ": Please specify a function for 'particleFactor'.");
  }
  if (m_2ndparticleFactor.isNull()) {
    throw gError("BondedPairParticleTensorNoiseVector::setup", "For " + m_properties.className() + ": Please specify a function for 'particleFactor'.");
  }
  if (m_pairFactor.isNull()) {
    throw gError("BondedPairParticleTensorNoiseVector::setup", "For " + m_properties.className() + ": Please specify a function for 'pairFactor'.");
  }
 
  m_pairFactor.setReturnType(Variant::VECTOR);
  m_pairFactor.addVariable("{dW}");
  m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_1stparticleFactor.addVariable("{dW}");
  m_2ndparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.addVariable("{dW}");

  if (m_randomize)
    m_rng.setSeed(getpid());
  else
    m_rng.setSeed(RNG_DEFAULT_SEED);

  BondedPairParticleCalc::setup();

  MSG_DEBUG("BondedPairParticleTensorNoiseVector::setup", ": END: m_symmetry = " << m_symmetry << ", &m_symmetry=" << &m_symmetry);    

 
}

void BondedPairParticleTensorNoiseVector::copyMembersTo(ValCalculator* vc)
{
  MSG_DEBUG("BondedPairParticleTensorNoiseVector::copyMembersTo", "CALLED!");

  ((BondedPairParticleTensorNoiseVector*) vc)->m_controller = m_controller;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_randomize = m_randomize;
//   ((BondedPairParticleTensorNoiseVector*) vc)->m_noise = m_noise;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_symmetry = m_symmetry;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_pairFactor = m_pairFactor;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_1stparticleFactor = m_1stparticleFactor;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_2ndparticleFactor = m_2ndparticleFactor;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_pairFactorStr = m_pairFactorStr;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_1stPExpression = m_1stPExpression;
  ((BondedPairParticleTensorNoiseVector*) vc)->m_2ndPExpression = m_2ndPExpression;

  BondedPairParticleCalc::copyMembersTo(vc);
}


#ifdef _OPENMP
void BondedPairParticleTensorNoiseVector::mergeCopies(ColourPair* cp, int thread_no) {

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
