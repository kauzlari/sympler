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



#include "f_pair_arbitrary.h"

#include "threads.h"
#include "particle.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

// const GenFTypeConcr<FPairArbitrary> f_pair_arbitrary("FPairArbitrary");

//---- Constructors/Destructor ----

FPairArbitrary::FPairArbitrary(Simulation *simulation): FPair(simulation)
{
  init();
}


FPairArbitrary::~FPairArbitrary()
{
}


void FPairArbitrary::init()
{
  m_properties.setClassName("FPairArbitrary");

  m_properties.setDescription( 
      "This is a completely general pair force FPV on a symbol x_i such that:\n"
      "dx_i = FPV*dt"
      "     = particleFactor_i*Sum_j(pairFactor_ij*weight_ij)*dt"
      "where pairFactor_ij includes all pair contributions of the pair ij, "
      "weight_ij represents the derivative of the used weighting function. "
      "particleFactor_i(j) represents factors specific to particle i(j).");

  STRINGPC  
      (pairFactor, m_pairFactorStr,
       "Function for pairFactor_ij. Type some nonsense to \n"
       "obtain a complete list of possible variables and constants.\n"
       "The expression may contain vectors and tensors,"
       " but as a whole it must represent the type of the corresponding force.");
  
  STRINGPC
      (particleFactor_i, m_1stPExpression,
       "The mathematical expression of the additional particle factor for the first particle."
      );

  STRINGPC
      (particleFactor_j, m_2ndPExpression,
       "The mathematical expression of the additional particle factor for the second particle."
      );

  INTPC
      (symmetry, m_symmetry, -2, 
       "'1' makes the interaction symmetric, '-1' antisymmetric."
      );
    
  m_symmetry = -1;

/*  m_pairFactorStr = "idVec(1)";
  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";*/
}

void FPairArbitrary::setup()
{
  FPair::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  
  assert(m_cp);
  
/*  if (m_pairFactor.isNull())
  throw gError("FPairVector::setup", "Please specify a function for 'pairFactor'.");*/
  m_pairFactor.setExpression(m_pairFactorStr);
//   m_pairFactor.setReturnType(Variant::VECTOR);
  m_pairFactor.setColourPair(m_cp);
  m_1stparticleFactor.setExpression(m_1stPExpression);
//   m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_1stparticleFactor.setColourPair(m_cp);
  m_2ndparticleFactor.setExpression(m_2ndPExpression);
//   m_2ndparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.setColourPair(m_cp);

  if (m_1stparticleFactor.isNull()) {
    throw gError("FPairArbitrary::setup", "For " + m_properties.className() + ": Please specify a function for 'particleFactor'.");
  }
  if (m_2ndparticleFactor.isNull()) {
    throw gError("FPairArbitrary::setup", "For " + m_properties.className() + ": Please specify a function for 'particleFactor'.");
  }
  if (m_pairFactor.isNull()) {
    throw gError("FPairArbitrary::setup", "For " + m_properties.className() + ": Please specify a function for 'pairFactor'.");
  }
  if(m_symmetry != 1 && m_symmetry != -1) 
    throw gError("FPairArbitrary::setup", "For " + m_properties.className() + ": m_symmetry = \"" + ObjToString(m_symmetry) + 
        "\" not allowed. Allowed values are '-1' and '1'.");  
}
