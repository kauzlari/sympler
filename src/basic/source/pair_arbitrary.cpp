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

#include "pair_arbitrary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


PairArbitrary::PairArbitrary(string symbol) 
  : ValCalculatorPair(symbol) {
}

PairArbitrary::PairArbitrary(Simulation* parent) 
  : ValCalculatorPair(parent) {
  init();
}

PairArbitrary::~PairArbitrary() {
}

void PairArbitrary::init() {
  
  STRINGPC
    (symbol, m_symbolName,
     "Name of the symbol.");
  
  STRINGPC(expression, m_expression_string,
	   "The mathematical expression to be computed.");

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Specifies the range of particle distances.");

  STRINGPC
    (useOldFor, m_oldSymbols,
     "Here, you can list the used symbols, which should be treated as \"old\", i.e., this calculator will not wait for those symbols to be computed beforehand, but it will take what it finds. Separate the symbols by the \"|\"- (\"pipe\"-) symbol."
     );

  m_expression_string == "undefined";
  
  m_oldSymbols = "---";
  
  m_cutoff = 0;

}

//---- Methods ----

void PairArbitrary::setup() {
  
  ValCalculatorPair::setup();
  
  if (m_cutoff <= 0)
    throw gError("PairArbitrary::setup for " + className(), ": Attribute 'cutoff' has no value > 0!");

  if (m_expression_string == "undefined")
    throw gError("PairArbitrary::setup for " + className(), ": Attribute 'expression' has value \"undefined\"");
  
  MSG_DEBUG("PairArbitrary::setup for " + className(), "  Cutoff used = " << m_cutoff);
  MSG_DEBUG("PairArbitrary::setup for " + className(), "  All pairs used = " << m_allPairs);
  MSG_DEBUG("PairArbitrary::setup for " + className(), "  Symbol used = " << m_symbolName);
  if (m_allPairs) {
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER,
       cp->setCutoff(m_cutoff);
       cp->setNeedPairs(true);
       MSG_DEBUG("PairArbitrary::setup for " + className() ," defined Cutoff = " << m_cutoff);
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
	throw gError("PairArbitrary::setup for " + className(), "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");
      if (m_species.second == "undefined")
	throw gError("PairArbitrary::setup for " + className(), "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");
      ColourPair* cp= M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);
      cp->setCutoff(m_cutoff);
      cp->setNeedPairs(true);
      m_function.setExpression(m_expression_string);
      m_function.setColourPair(cp);
      
    }
}


void PairArbitrary::setSlot(ColourPair* cp, size_t& slot,
						 bool oneProp) {  
  m_slot = slot = cp->tagFormat().addAttribute
    (m_symbolName + cp->toString(), m_datatype, false, m_symbolName).offset;
}

bool PairArbitrary::findStage()
{
  return Symbol::findStageNewPrelim();
}


bool PairArbitrary::findStage_0()
{
  return Symbol::findStageNewPrelim_0();
}


void PairArbitrary::addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
{
  FunctionParser::addToTypedValueList(m_function.usedSymbols(), usedSymbols);
}
