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


#include "bonded_pair_arbitrary.h"
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


BondedPairArbitrary::BondedPairArbitrary(Simulation* parent)
  : ValCalculatorPair(parent)
{
  init();
}

BondedPairArbitrary::~BondedPairArbitrary()
{
}

void BondedPairArbitrary::init()
{
//   m_properties.setClassName("BondedPairArbitrary");
  m_properties.setClassName("ValCalculatorPair");
  // Note that if a child-class uses setClassName, it also sets back the name. Since this is the usual case, the next line is rather a reminder that the children have to do this line again according to their needs!
  m_properties.setName("BondedPairArbitrary");

  STRINGPC(listName, m_listName, "Identifier of the list of bonded pairs, this Calculator should belong to.");

  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the calculated property.");

  STRINGPC
      (expression, m_expression,
       "The mathematical expression to be computed for the pairFactor."
      );

  STRINGPC
      (useOldFor, m_oldSymbols,
       "Here, you can list the used symbols, which should be treated as \"old\", i.e., this calculator will not wait for those symbols to be computed beforehand, but it will take what it finds. Separate the symbols by the \"|\"- (\"pipe\"-) symbol."
      );

#ifdef _OPENMP
  m_particleCalculator = false;
#endif

  m_listName = "undefined";

  m_oldSymbols = "---";

}

void BondedPairArbitrary::setup()
{
   if(m_allPairs == true)
     throw gError("BondedPairArbitrary::setup", "For module " + className() + " for symbol " + m_symbolName + ": For bonded pairs, 'allPairs = \"yes\" is currently not supported. Please switch off.");

  if(m_listName == "undefined")
    throw gError("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": Attribute 'listName' is undefined!");

  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": Attribute 'stage' is " + ObjToString(m_phaseUser) + ", which is none of the allowed values \"0\", \"1\", \"2\".");

  /*m_allPairs == false ALWAYS!*/
   
  if(m_species.first == "undefined")
    throw gError("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": Attribute 'species1' has value \"undefined\"."); 
  if(m_species.second == "undefined")
    throw gError("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": Attribute 'species2' has value \"undefined\"."); 
  
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);
    
  // next lines are just necessary for non-bonded pairs, so not here!
  //      cp->setNeedPairs(true);
  //   cp->setCutoff(m_cutoff);

  if(m_expression == "undefined")
    throw gError("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": Attribute 'expression' has value \"undefined\"");

  m_function.setExpression(m_expression);
  m_function.setColourPair(cp);

  
  MSG_DEBUG("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": registering " << this->name() << ".");
  
  if(m_phaseUser == 0)    
    cp->registerBondedCalc_0(this);
  else if(m_phaseUser == 1)    
    cp->registerBondedCalc(this);
  else // so it is 2
    {
      ValCalculator* vc = copyMySelf() /*new BondedPairArbitrary(*this)*/;

      copyMembersTo(vc);
      
      MSG_DEBUG("BondedPairArbitrary::setup", "Module " + className() + " for symbol " + m_symbolName + ": registering copy of " << this->name() << ", CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");    
      
      cp->registerBondedCalc(vc);
      cp->registerBondedCalc_0(this);
    }

  ValCalculatorPair::setup();

}

void BondedPairArbitrary::setSlot(ColourPair* cp, size_t& slot,
						 bool oneProp) {
  if(!cp->tagFormat().attrExists(m_symbolName))
    m_slot = slot = cp->tagFormat().addAttribute
      (m_symbolName, m_datatype, false, m_symbolName).offset;
  else
    throw gError("BondedPairArbitrary::setSlot for " + className(), "Symbol " + m_symbolName + " already exists for species pair (" + cp->firstSpecies() + ", " + cp->secondSpecies() + "). This would mean this module overwrites any previous values. So we should not do this.");
}

void BondedPairArbitrary::copyMembersTo(ValCalculator* vc)
{

  ((BondedPairArbitrary*) vc)->m_overwrite = m_overwrite;

#ifdef _OPENMP
  ((BondedPairArbitrary*) vc)->m_particleCalculator = m_particleCalculator;
#endif

   ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  ((BondedPairArbitrary*) vc)->m_function.setExpression(m_expression);
  ((BondedPairArbitrary*) vc)->m_function.setColourPair(cp);

}


void BondedPairArbitrary::addMyUsedSymbolsTo(typed_value_list_t& usedSymbols) {

  FunctionParser::addToTypedValueList(m_function.usedSymbols(), usedSymbols);      
}


