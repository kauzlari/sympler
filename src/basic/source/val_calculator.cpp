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



#include "manager_cell.h"
#include "val_calculator.h"
#include "simulation.h"
#include "colour_pair.h"
/*#include "val_calculator_rho.h"
#include "val_calculator_volume.h"*/
// #include "val_calculator_kernel.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculator::ValCalculator(string symbol): Symbol(symbol)
{
//   m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  //   MSG_DEBUG("ValCalculator::ValCalculator", "CONSTRUCTOR");
}

ValCalculator::ValCalculator(/*Node*/Simulation* parent): Symbol(parent)
{
//   m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  init();
//   MSG_DEBUG("ValCalculator::ValCalculator", "NODE-CONSTRUCTOR");
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
  MSG_DEBUG("ValCalculator::cleanPairSymbol", className() << ": shortened name of symbol: " << name);
}

ValCalculatorPair::ValCalculatorPair(/*Node*/Simulation* parent): ValCalculator(parent)/*, m_persistency(false)*/
{
  init();
//   MSG_DEBUG("ValCalculator::ValCalculator", "CONSTRUCTOR");
}

void ValCalculatorPair::init()
{
  m_properties.setClassName("ValCalculatorPair");

  BOOLPC(allPairs, m_allPairs, "Should the quantity be computed for all colour combinations? If yes, this disables 'species1' and 'species2'. NOTE: This feature is currently not supported for bonded pairs.");

  m_allPairs = false;

}

void ValCalculatorPair::setup()
{
  
  if(m_allPairs)
  {
    FOR_EACH_COLOUR_PAIR
    (
      M_MANAGER,
      if(cp->tagFormat().attrExists(m_symbolName))
        throw gError("ValCalculatorPair::setup", "Symbol " + m_symbolName + " already existing. Second definition is not allowed for this Calculator.");

      // see CONVENTION5 for rule about persistencies
#if 0
      // new rules
      // FIXME: not checked how meaningful this is for the case that both species don't have an IntegratorPosition
      // by the way, how meaningful is this case itself ?!?
      if (!M_CONTROLLER->findIntegrator("IntegratorPosition", cp->firstSpecies()) && 
           !M_CONTROLLER->findIntegrator("IntegratorPosition", cp->secondSpecies()))
        m_persistency = true;
#endif
      m_slot = cp->tagFormat().addAttribute(m_symbolName, m_datatype, false/*persistency*//*m_persistency*/, m_symbolName).offset;

      vector<ColourPair*>::iterator cpTester = __cp;
    // is it the last calculator to be created?
      if(++cpTester == __end)
        cp -> registerCalc(this);
        // No? Then make a copy
      else 
      {
        ValCalculator* vc = copyMySelf()/*new ValCalculatorPair(*this)*/;
        assert(((ValCalculatorPair*) vc)->m_symbolName == m_symbolName);
        assert(((ValCalculatorPair*) vc)->stage() == m_stage);
//         assert(vc->m_wf == m_wf);
//         assert(vc->m_wfName == m_wfName);
        assert(((ValCalculatorPair*) vc)->m_slot == m_slot);
        assert(((ValCalculatorPair*) vc)->m_parent == m_parent);
        assert(((ValCalculatorPair*) vc)->m_species.first == m_species.first);
        assert(((ValCalculatorPair*) vc)->m_species.second == m_species.second);

        cp->registerCalc(vc);
      }
    );    
  }
  else
  {
    if(m_species.first == "undefined")
      throw gError("ValCalculatorPair::setup", "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled."); 
    if(m_species.second == "undefined")
      throw gError("ValCalculatorPair::setup", "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled."); 

    ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

    // see CONVENTION5 for rule about persistencies
#if 0
    // new rules
    // FIXME: not checked how meaningful this is for the case that both species don't have an IntegratorPosition
    // by the way, how meaningful is this case itself ?!?
    if (!M_CONTROLLER->findIntegrator("IntegratorPosition", cp->firstSpecies()) && 
         !M_CONTROLLER->findIntegrator("IntegratorPosition", cp->secondSpecies()))
      m_persistency = true;
#endif
      
    if(cp->tagFormat().attrExists(m_symbolName))
      throw gError("ValCalculatorPair::setup", "Symbol " + m_symbolName + " already existing. Second definition is not allowed for this Calculator.");
      
    m_slot = cp->tagFormat().addAttribute(m_symbolName, m_datatype, false/*persistency*//*m_persistency*/, m_symbolName).offset;
    
    cp -> registerCalc(this);
  }

}

void /*pair<size_t, size_t>*/ ValCalculatorPair::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
{
  m_slot = slot = cp->tagFormat().addAttribute
      ("ValCalculator_" + myName() + "_" + cp->toString(), m_datatype).offset;
}

ValCalculatorPart::ValCalculatorPart(/*Node*/Simulation* parent): ValCalculator(parent)
{
  init();
//   MSG_DEBUG("ValCalculator::ValCalculator", "CONSTRUCTOR");
}

void ValCalculatorPart::init()
{
  m_properties.setClassName("ValCalculatorPart");
 
#ifdef _OPENMP
  m_copy_slots.resize(global::n_threads);
#endif

}

void ValCalculatorPart::setup()
{

  throw gError("ValCalculatorPart::setup", "Shouldn't currently be called. Contact the programmer.");
  
}


NonBondedPairParticleCalculator::NonBondedPairParticleCalculator(string symbol)
  : ValCalculatorPart(symbol)
{

}

NonBondedPairParticleCalculator::NonBondedPairParticleCalculator(/*Node*/Simulation* parent): ValCalculatorPart(parent)
{
  init();
}


void NonBondedPairParticleCalculator::init()
{
  BOOLPC(allPairs, m_allPairs, "Should the quantity be computed for all colour combinations? If yes, this disables 'species1' and 'species2'");

  m_allPairs = false;


}


