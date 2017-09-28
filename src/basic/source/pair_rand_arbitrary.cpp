/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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



#include "pair_rand_arbitrary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


PairRandArbitrary::PairRandArbitrary(string symbol) :
	ValCalculatorPair(symbol) {
}

PairRandArbitrary::PairRandArbitrary(Simulation* parent) :
	ValCalculatorPair(parent) {
	init();
}

PairRandArbitrary::~PairRandArbitrary() {
}

void PairRandArbitrary::init() {
  m_properties.setClassName("PairRandArbitrary");
  
  STRINGPC
    (symbol, m_symbolName,
     "Name of the random symbol.")
    ;
  
  STRINGPC(pairFactor, m_expression_string,
	   "The mathematical expression to be computed for the pairFactor. It may contain scalars, vectors and matrices, but in the end it must ve of the same type as the final symbol. In the case of vectors or matrices, component-wise multiplication is performed.")
    ;

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Specifies the maximum distance of a pair of particles. Beyond this distance, this module leaves the to be computed symbol zero.")
    ;
  
  STRINGPC
    (useOldFor, m_oldSymbols,
     "Here, you can list the used symbols in the pair-factor, which should be treated as \"old\", i.e., this calculator will not wait for those symbols to be computed beforehand, but it will take what it finds. Separate the symbols by the \"|\"- (\"pipe\"-) symbol."
     )
    ;
  
  INTPC
    (seed, m_seed, 0,
     "Seed to be used for the random number generator."
     )
    ;
  
  m_oldSymbols = "---";
  
  m_cutoff = 0;  

  m_seed = RNG_DEFAULT_SEED;
  
}

//---- Methods ----

void PairRandArbitrary::setup() {
  
  ValCalculatorPair::setup();

  m_rng.setSeed(m_seed);
  
  if (m_cutoff <= 0)
    throw gError("PairRandArbitrary::setup", className() + ": Attribute 'cutoff' has no value > 0!");

  if (m_expression_string == "undefined")
    throw gError("PairRandArbitrary::setup", className() + ": Attribute 'pairFactor' has value \"undefined\"");
  
  MSG_DEBUG("PairRandArbitrary::setup",className() << "  Cutoff used  = " << m_cutoff);

  if (m_allPairs) {
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER,
       cp->setCutoff(m_cutoff);
       cp->setNeedPairs(true);
       vector<ColourPair*>::iterator cpTester = __cp;
       if(++cpTester == __end)
	 {
	   m_function.setExpression(m_expression_string);
	   m_function.setColourPair(cp);
	 }
       )
      ;
  } else /*m_allPairs == false*/
    {
      if (m_species.first == "undefined")
	throw gError("PairRandArbitrary::setup", "For module " + className() + ": Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");
      if (m_species.second == "undefined")
	throw gError("PairRandArbitrary::setup", "For module " + className() + ": Attribute 'species2' has value \"undefined\" and 'allPairs' is disabled.");
      ColourPair* cp= M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);
      cp->setCutoff(m_cutoff);
      cp->setNeedPairs(true);
      m_function.setExpression(m_expression_string);
      m_function.setColourPair(cp);
      
    }
  
  
}

void PairRandArbitrary::setSlot(ColourPair* cp, size_t& slot, bool oneProp) {
  m_slot = slot = cp->tagFormat().addAttribute
    (m_symbolName, m_datatype, false, m_symbolName).offset;
}


bool PairRandArbitrary::findStage()
{
  return Symbol::findStageNewPrelim();
}


bool PairRandArbitrary::findStage_0()
{
  return Symbol::findStageNewPrelim_0();
}


void PairRandArbitrary::addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
{
  FunctionParser::addToTypedValueList(m_function.usedSymbols(), usedSymbols);
}
