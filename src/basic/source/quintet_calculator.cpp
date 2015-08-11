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


#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "quintet_calculator.h"

#include "quintet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

//---- Constructors/Destructor ----


QuintetCalculator::QuintetCalculator(Simulation *simulation): Symbol(simulation)
{
  // does not depend on other symbols
  m_stage = 0;

  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void QuintetCalculator::init()
{
  m_properties.setClassName("QuintetCalcPart");
  m_properties.setName("QuintetCalculator");

  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the calculated property.");
  
  STRINGPC(listName, m_listName, "Name of the quintet list to compute on.");

  STRINGPC(species1, m_species[0], "Species of first particle in the quintets.");
  STRINGPC(species2, m_species[1], "Species of second particle in the quintets.");
  STRINGPC(species3, m_species[2], "Species of third particle in the quintets.");
  STRINGPC(species4, m_species[3], "Species of second particle in the quintets.");
  STRINGPC(species5, m_species[4], "Species of third particle in the quintets.");
  
  BOOLPC(periodic, m_periodic, "Should periodic boundary conditions be applied to the connections?")

  BOOLPC
      (overwrite, m_overwrite,
       "Is this calculator allowed to overwrite already existing symbols " 
           "with name 'symbol' ?");

  m_symbolName = "undefined";

  m_listName = "undefined";

  m_species[0] = "undefined";
  m_species[1] = "undefined";
  m_species[2] = "undefined";
  m_species[3] = "undefined";
  m_species[4] = "undefined";

  m_overwrite = false;
  m_periodic = true;
}


void QuintetCalculator::setup()
{
  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

  if(m_symbolName == "undefined")
    throw gError("QuintetCalculator::setup", "Attribute 'symbol' undefined!");

  if(m_listName == "undefined")
    throw gError("QuintetCalculator::setup", "Attribute 'listName' undefined!");

  // add attributes


    if(m_species[0] == "undefined")
      throw gError("QuintetCalculator::setup", ": Attribute 'species1' has value \"undefined\"!"); 
    if(m_species[1] == "undefined")
      throw gError("QuintetCalculator::setup", ": Attribute 'species2' has value \"undefined\"!"); 
    if(m_species[2] == "undefined")
      throw gError("QuintetCalculator::setup", ": Attribute 'species3' has value \"undefined\"!"); 

    m_firstColour = M_MANAGER->getColour(m_species[0]);
    m_secondColour = M_MANAGER->getColour(m_species[1]);
    m_thirdColour = M_MANAGER->getColour(m_species[2]);
    m_fourthColour = M_MANAGER->getColour(m_species[3]);
    m_fifthColour = M_MANAGER->getColour(m_species[4]);


    if(m_overwrite)
    {

      for(size_t p = 0; p < 3; ++p) {
	try
	  {
	    m_slots[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].indexOf(m_symbolName, m_datatype);
	    m_slots[p] = Particle::s_tag_format[M_MANAGER->getColour(m_species[p])].offsetByIndex(m_slots[p]);
	  }
	catch(gError& err)
	  {
	    throw gError("QuintetCalculator::setup", ": search for symbol for species '" + m_species[p] + " failed. The message was " + err.message()); 
	  }
      } // end for(size_t p = 0

    }
    else // m_overwrite = false
    {
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].attrExists(m_symbolName)) 
        throw gError("QuintetCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[0] + "'. Second definition is not allowed for overwrite = \"no\"");

      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].attrExists(m_symbolName)) 
        throw gError("QuintetCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[1] + "'. Second definition is not allowed for overwrite = \"no\"");

      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].attrExists(m_symbolName)) 
        throw gError("QuintetCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[2] + "'. Second definition is not allowed for overwrite = \"no\"");

   if(Particle::s_tag_format[M_MANAGER->getColour(m_species[3])].attrExists(m_symbolName)) 
        throw gError("QuintetCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[3] + "'. Second definition is not allowed for overwrite = \"no\"");

   if(Particle::s_tag_format[M_MANAGER->getColour(m_species[4])].attrExists(m_symbolName)) 
        throw gError("QuintetCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[4] + "'. Second definition is not allowed for overwrite = \"no\"");

       MSG_DEBUG("QuintetCalculator::setup", "adding symbol with name " << m_symbolName);
      
      // see CONVENTION5 for rule about persistencies
      m_slots[0] = Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
  
      if(m_species[0] != m_species[1])
      {
        // see CONVENTION5 for rule about persistencies
	m_slots[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
      }
      else m_slots[1] = m_slots[0];

      if(m_species[2] == m_species[0])
	m_slots[2] = m_slots[0];
      else if(m_species[2] == m_species[1])
	m_slots[2] = m_slots[1];
      else
        // see CONVENTION5 for rule about persistencies
	m_slots[2] = Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;

    } // end of else of if(m_overwrite == false)        
    
    if(m_phaseUser == 0)    
      M_PHASE->registerBondedQuinCalc_0(this);
    else if(m_phaseUser == 1)    
      M_PHASE->registerBondedQuinCalc(this);
    else // so it is 2
    {
      QuintetCalculator* vc = copyMySelf();
      assert(((QuintetCalculator*) vc)->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((QuintetCalculator*) vc)->m_slots[0] == m_slots[0]);
      assert(((QuintetCalculator*) vc)->m_slots[1] == m_slots[1]);
      assert(((QuintetCalculator*) vc)->m_slots[2] == m_slots[2]);
      assert(((QuintetCalculator*) vc)->m_slots[3] == m_slots[3]);
      assert(((QuintetCalculator*) vc)->m_slots[4] == m_slots[4]);

      assert(((QuintetCalculator*) vc)->m_parent == m_parent);
      assert(((QuintetCalculator*) vc)->m_overwrite == m_overwrite);
      assert(((QuintetCalculator*) vc)->m_species[0] == m_species[0]);
      assert(((QuintetCalculator*) vc)->m_species[1] == m_species[1]);
      assert(((QuintetCalculator*) vc)->m_species[2] == m_species[2]);
      assert(((QuintetCalculator*) vc)->m_species[3] == m_species[3]);
      assert(((QuintetCalculator*) vc)->m_species[4] == m_species[4]);
        
      MSG_DEBUG("QuintetCalculator::setup", ": registering copy.");    
       
      M_PHASE->registerBondedQuinCalc(vc);
      M_PHASE->registerBondedQuinCalc_0(this);
    }
}

void QuintetCalculator::setupAfterParticleCreation()
{
  m_boxSize = M_PHASE->boundary()->boundingBox().size();

  // FOLLOWING check cannot be done here because quintet lists 
  // are also created in a setupAfterParticleCreation()

//   // check for consistency of calculator's and quintets' species

//   quintetList* trl = M_PHASE->returnQuintetList(m_listName);

//   quintetList::iterator firstTr = trl->begin();

//   if(M_MANAGER->getColour(m_species[0]) != firstTr->a->c)
//     throw gError("QuintetCalculator::setupAfterParticleCreation", "Inconsistency between the attribute 'species1' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[1]) != firstTr->b->c)
//     throw gError("QuintetCalculator::setupAfterParticleCreation", "Inconsistency between the attribute 'species2' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[2]) != firstTr->c->c)
//     throw gError("QuintetCalculator::setupAfterParticleCreation", "Inconsistency between the attribute'species3' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");

}
/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool QuintetCalculator::findStage()
{
  
  // currently (2009/07/28) this always returns true
  return Symbol::findStage();
}


// check for consistency of calculator's and quintets' species
void QuintetCalculator::checkConsistency() {

  quintetList* quin = M_PHASE->returnQuintetList(m_listName);
  
  quintetList::iterator firstQuin = quin->begin();
  
  if(M_MANAGER->getColour(m_species[0]) != firstQuin->p00->c)
    throw gError("QuintetCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species1=\"" + m_species[0] + "\" of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  if(M_MANAGER->getColour(m_species[1]) != firstQuin->p02->c)
    throw gError("QuintetCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species2=\"" + m_species[1] + "\"  of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  if(M_MANAGER->getColour(m_species[2]) != firstQuin->p22->c)
    throw gError("QuintetCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species3=\"" + m_species[2] + "\"  of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  if(M_MANAGER->getColour(m_species[3]) != firstQuin->p20->c)
    throw gError("QuintetCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species4=\"" + m_species[3] + "\"  of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  if(M_MANAGER->getColour(m_species[4]) != firstQuin->p11->c)
    throw gError("QuintetCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species5=\"" + m_species[4] + "\"  of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");

}   


/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool QuintetCalculator::findStage_0()
{
//   // Following check done here because only suitable function after setupAfterParticleCreation() (see also comment there)
//   // check for consistency of calculator's and quintets' species
  
//   quintetList* trl = M_PHASE->returnQuintetList(m_listName);
  
//   quintetList::iterator firstTr = trl->begin();
  
//   if(M_MANAGER->getColour(m_species[0]) != firstTr->a->c)
//     throw gError("QuintetCalculator::findStage_0", "Inconsistency between the attribute 'species1' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[1]) != firstTr->b->c)
//     throw gError("QuintetCalculator::findStage_0", "Inconsistency between the attribute 'species2' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[2]) != firstTr->c->c)
//     throw gError("QuintetCalculator::findStage_0", "Inconsistency between the attribute'species3' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  
  // currently (2009/07/28) this always returns true
  return Symbol::findStage_0();
}
