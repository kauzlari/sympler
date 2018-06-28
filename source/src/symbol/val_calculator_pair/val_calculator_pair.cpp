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


#include "manager_cell.h"
#include "val_calculator.h"
#include "val_calculator_pair.h"
#include "simulation.h"
#include "colour_pair.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;


ValCalculatorPair::ValCalculatorPair(/*Node*/Simulation* parent)
  : ValCalculator(parent)
{
  init();
}

void ValCalculatorPair::init()
{
  m_properties.setClassName("ValCalculatorPair");

  BOOLPC(allPairs, m_allPairs, "Should the quantity be computed for "
	 "all colour combinations? If yes, this disables 'species1' "
	 "and 'species2'. NOTE: This feature is currently not "
	 "supported for bonded pairs.");

  m_allPairs = false;

}

void ValCalculatorPair::setup()
{
  
  if(m_allPairs) {
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER,
       if(cp->tagFormat().attrExists(m_symbolName))
	 throw
	   gError("ValCalculatorPair::setup", "Symbol " + m_symbolName
		  + " already existing. Second definition is not "
		  "allowed for this Calculator.");
       
       // see CONVENTION5 for rule about persistencies
       m_slot =
       cp->tagFormat().addAttribute
       (m_symbolName, m_datatype, false/*persistency*/, m_symbolName)
       .offset;
       
       vector<ColourPair*>::iterator cpTester = __cp;
       
       // is it the last calculator to be created?
       if(++cpTester == __end) {
	 if(m_phaseUser == 0)
	   cp -> registerCalc_0(this);
	 else if(m_phaseUser == 1)
	   cp -> registerCalc(this);
	 // so it is 2
	 else {
	   ValCalculator* vc = copyMySelf();
	   makeCopySafe((ValCalculatorPair*) vc);
	   cp -> registerCalc_0(vc);
	   cp -> registerCalc(this);
	 }
       } // end of if(++cpTester == __end)
       // No? Then make a copy in any case
       else {
	 ValCalculator* vc = copyMySelf();
	 makeCopySafe((ValCalculatorPair*) vc);
	 
	 if(m_phaseUser == 0)
	   cp -> registerCalc_0(vc);
	 else if(m_phaseUser == 1)
	   cp -> registerCalc(vc);
	 // so it is 2
	 else {
	   ValCalculator* vc2 = copyMySelf();
	   makeCopySafe((ValCalculatorPair*) vc2);
	   cp -> registerCalc_0(vc);
	   cp -> registerCalc(vc2);
	 }
       } // end of else of if(++cpTester == __end)
       ); // end of FOR_EACH_COLOUR_PAIR
    
  } // end of if(m_allPairs)
  else {
    
    if(m_species.first == "undefined")
      throw
	gError("ValCalculatorPair::setup", "Attribute 'species1' has "
	       "value \"undefined\" and 'allPairs' is disabled."); 

    if(m_species.second == "undefined")
      throw
	gError("ValCalculatorPair::setup", "Attribute 'species1' has "
	       "value \"undefined\" and 'allPairs' is disabled."); 

    ColourPair* cp =
      M_MANAGER->cp
      (M_MANAGER->getColour(m_species.first),
       M_MANAGER->getColour(m_species.second));

    // see CONVENTION5 for rule about persistencies      
    if(cp->tagFormat().attrExists(m_symbolName))
      throw
	gError("ValCalculatorPair::setup", "Symbol " + m_symbolName +
	       " already existing. Second definition is not allowed "
	       "for this Calculator.");
      
    m_slot =
      cp->tagFormat().addAttribute
      (m_symbolName, m_datatype, false/*persistency*/, m_symbolName)
      .offset;

    if(m_phaseUser == 0)
      cp -> registerCalc_0(this);
    else if(m_phaseUser == 1)
      cp -> registerCalc(this);
    // so it is 2
    else {
      ValCalculator* vc = copyMySelf();
      makeCopySafe((ValCalculatorPair*) vc);
      cp -> registerCalc_0(vc);
      cp -> registerCalc(this);
    }

  } // end of else of if(m_allPairs)

}

void ValCalculatorPair::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
{
  m_slot = slot = cp->tagFormat().addAttribute
    ("ValCalculator_" + myName() + "_" + cp->toString(), m_datatype
     /*, persistency = false, symbol = ""*/).offset;
}

void ValCalculatorPair::makeCopySafe(ValCalculatorPair* vc) const {
  assert(((ValCalculatorPair*) vc)->m_symbolName == m_symbolName);
  assert(((ValCalculatorPair*) vc)->stage() == m_stage);
  assert(((ValCalculatorPair*) vc)->m_slot == m_slot);
  assert(((ValCalculatorPair*) vc)->m_parent == m_parent);
  assert(((ValCalculatorPair*) vc)->m_species.first == m_species.first);
  assert(((ValCalculatorPair*) vc)->m_species.second == m_species.second);
}


