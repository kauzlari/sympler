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
#include "simulation.h"
#include "colour_pair.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculator::ValCalculator(string symbol): Symbol(symbol)
{
  m_datatype = DataFormat::DOUBLE;
}

ValCalculator::ValCalculator(/*Node*/Simulation* parent): Symbol(parent)
{
  m_datatype = DataFormat::DOUBLE;
  init();
}

void ValCalculator::init()
{
  m_properties.setClassName("ValCalculator");

  STRINGPC
      (species1, m_species.first,
       "Name for the first species of the pairs, this Symbol is used for.");
  
  STRINGPC
      (species2, m_species.second,
       "Name for the second species of the pairs, this Symbol is used for.");
  
  m_species.first = "undefined";
  m_species.second = "undefined"; 
}

void ValCalculator::setup()
{
  Symbol::setup();
}

void ValCalculator::cleanSymbol(string& name) const
{
  if(name[0] == '{' || name[0] == '[') {
    // remove the first bracket
    name.erase(0, 1);
    // remove the last bracket; don't know why, but with these arguments it works
    name.erase(name.size()-1, name.size()-1);
  }
  // and the "i", "j" and "ij" of the pair expression have to be removed
  if(name[name.size()-2] == 'i' && name[name.size()-1] == 'j')
    // remove "ij"
    name.erase(name.size()-2, name.size()-1); 
  else
    // remove "i" or "j"
    name.erase(name.size()-1, name.size()-1);
  
}

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

ValCalculatorPart::ValCalculatorPart(/*Node*/Simulation* parent)
  : ValCalculator(parent)
{
  init();
}

void ValCalculatorPart::init()
{
  m_properties.setClassName("ValCalculatorPart");

  // START: unfinished stuff from 2018-05-08 ////////////////////
  
  // BOOLPC
  //   (selfReset, m_selfReset, "Should this module protect its computed "
  //    "symbol from automatic resetting (to zero) and reset it by "
  //    "itself? The self-reset will be done"
  //    "immediately before the start of computations by any Symbols, "
  //    "including this one in stage 0 and/or 1 as selected. This may be "
  //    "useful to prevent a too early reset, for example by a triggered "
  //    "position update or neighbour list rebuild. In selfReset "
  //    "mode, overwriting existing symbols "
  //    "('overrite = \"yes\"') by this module is forbidden.");

  // m_selfReset = false;

  // END: unfinished stuff from 2018-05-08 ////////////////////

  
#ifdef _OPENMP
  m_copy_slots.resize(global::n_threads);
#endif

}

void ValCalculatorPart::setup()
{
  ValCalculator::setup();

  // START: unfinished stuff from 2018-05-08 ////////////////////
  
  // if(m_selfReset) {
  //   if(m_overwrite)
  //     // Not all children allow to control m_overwrite via XML-input. 
  //     // But since the default in class Symbol is m_overwrite = false,
  //     // this exception and its message should always make sense when
  //     // thrown, since m_overwrite = true is either set via XML or via
  //     // an error in the implementation of a child class. The latter
  //     // should only be metioned in this comment and *not* in the
  //     // exception message, such that we don't confuse users.
  //     throw gError("ValCalculatorPart::setup for module " + className(),
  // 		   "'overwrite = \"yes\"' is not allowed together with "
  // 		   "'selfReset = \"yes\"'.");
  //   m_persistency = true;

  //   // The resetting will be done by an overriden Node::precompute()
  //   if(m_phase == 0) M_CONTROLLER->registerForPrecomputation_0(this);
  //   if(m_phase == 1) M_CONTROLLER->registerForPrecomputation(this);
  //   if(m_phase == 2) {
  //     M_CONTROLLER->registerForPrecomputation(this);
  //     M_CONTROLLER->registerForPrecomputation_0(this);
  //   }
  // } // end of if if(m_selfReset)
  // else
  //   m_persistency = false;

  // END: unfinished stuff from 2018-05-08 ////////////////////
  
}


NonBondedPairParticleCalculator::NonBondedPairParticleCalculator(string symbol)
  : ValCalculatorPart(symbol)
{

}


NonBondedPairParticleCalculator::NonBondedPairParticleCalculator
(/*Node*/Simulation* parent)
  : ValCalculatorPart(parent)
{
  init();
}


void NonBondedPairParticleCalculator::init()
{
  BOOLPC(allPairs, m_allPairs,
	 "Should the quantity be computed for all colour combinations? "
	 "If \"yes\", this disables 'species1' and 'species2'");
  
  m_allPairs = false;
}


