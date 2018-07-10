/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


#include "bonded_pair_particle_arbitrary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"



#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


BondedPairParticleArbitrary::BondedPairParticleArbitrary(Simulation* parent)
  : BondedPairParticleCalc(parent)
{
  init();
}

BondedPairParticleArbitrary::~BondedPairParticleArbitrary()
{
}

void BondedPairParticleArbitrary::init()
{
//   m_properties.setClassName("BondedPairParticleArbitrary");
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("BondedPairParticleArbitrary");

  STRINGPC
      (expression, m_expression,
       "The mathematical expression to be computed for the pairFactor."
      );
  
  STRINGPC
      (particleFactor_i, m_1stPExpression,
       "The mathematical expression of the additional particle factor for the first particle."
      );

  STRINGPC
      (particleFactor_j, m_2ndPExpression,
       "The mathematical expression of the additional particle factor for the second particle."
      );

  STRINGPC
      (useOldFor, m_oldSymbols,
       "Here, you can list the used symbols, which should be treated as \"old\", i.e., this calculator will not wait for those symbols to be computed beforehand, but it will take what it finds. Separate the symbols by the \"|\"- (\"pipe\"-) symbol."
      );
    
  m_oldSymbols = "---";

}

void BondedPairParticleArbitrary::setup()
{
  if(m_expression == "undefined")
    throw gError("BondedPairParticleArbitrary::setup", ": Attribute 'expression' has value \"undefined\"");
       
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  
  m_function.setExpression(m_expression);
  m_function.setColourPair(cp);
  m_1stparticleFactor.setExpression(m_1stPExpression);
  m_1stparticleFactor.setColourPair(cp);
  m_2ndparticleFactor.setExpression(m_2ndPExpression);
  m_2ndparticleFactor.setColourPair(cp);

  BondedPairParticleCalc::setup();
  
}


void BondedPairParticleArbitrary::copyMembersTo(ValCalculator* vc)
{
  BondedPairParticleCalc::copyMembersTo(vc);

   ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  ((BondedPairParticleArbitrary*) vc)->m_function.setExpression(m_expression);
  ((BondedPairParticleArbitrary*) vc)->m_function.setColourPair(cp);
  //       ((BondedPairParticleArbitrary*) vc)->m_function.setReturnType(funcType);
  ((BondedPairParticleArbitrary*) vc)->m_1stparticleFactor.setExpression(m_1stPExpression);
  ((BondedPairParticleArbitrary*) vc)->m_1stparticleFactor.setColourPair(cp);
  ((BondedPairParticleArbitrary*) vc)->m_2ndparticleFactor.setExpression(m_2ndPExpression);
  ((BondedPairParticleArbitrary*) vc)->m_2ndparticleFactor.setColourPair(cp);

}


void BondedPairParticleArbitrary::addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
{
  FunctionParser::addToTypedValueList(m_function.usedSymbols(), usedSymbols);
  FunctionParser::addToTypedValueList(m_1stparticleFactor.usedSymbols(), usedSymbols);
  FunctionParser::addToTypedValueList(/*const typed_value_list_t&*/ m_2ndparticleFactor.usedSymbols(), /*typed_value_list_t&*/ usedSymbols);
}

