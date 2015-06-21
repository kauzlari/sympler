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


#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "triplet_calculator.h"

#include "triplet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

//---- Constructors/Destructor ----


TripletCalculator::TripletCalculator(Simulation *simulation): Symbol(simulation)
{
  // does not depend on other symbols
  m_stage = 0;

  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void TripletCalculator::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalculator");

  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the calculated property.");
  
  STRINGPC(listName, m_listName, "Name of the triplet list to compute on.");

  STRINGPC(species1, m_species[0], "Species of first particle in the triplets.");
  STRINGPC(species2, m_species[1], "Species of second particle in the triplets.");
  STRINGPC(species3, m_species[2], "Species of third particle in the triplets.");
  
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

  m_overwrite = false;
  m_periodic = true;
}


void TripletCalculator::setupTags()
{
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
	    throw gError("TripletCalculator::setup", ": search for symbol for species '" + m_species[p] + " failed. The message was " + err.message()); 
	  }
      } // end for(size_t p = 0
      
    }
  else // m_overwrite = false
    {
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].attrExists(m_symbolName)) 
        throw gError("TripletCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[0] + "'. Second definition is not allowed for overwrite = \"no\"");
      
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].attrExists(m_symbolName)) 
        throw gError("TripletCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[1] + "'. Second definition is not allowed for overwrite = \"no\"");
      
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].attrExists(m_symbolName)) 
        throw gError("TripletCalculator::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[2] + "'. Second definition is not allowed for overwrite = \"no\"");
      
      
      MSG_DEBUG("TripletCalculator::setup", "adding symbol with name " << m_symbolName);
      
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

}


void TripletCalculator::registerCalc()
{
  if(m_phaseUser == 0)    
    M_PHASE->registerBondedCalc_0(this);
  else if(m_phaseUser == 1)    
    M_PHASE->registerBondedCalc(this);
  else // so it is 2
    {
      TripletCalculator* vc = copyMySelf();
      assert(((TripletCalculator*) vc)->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((TripletCalculator*) vc)->m_slots[0] == m_slots[0]);
      assert(((TripletCalculator*) vc)->m_slots[1] == m_slots[1]);
      assert(((TripletCalculator*) vc)->m_slots[2] == m_slots[2]);
      assert(((TripletCalculator*) vc)->m_parent == m_parent);
      assert(((TripletCalculator*) vc)->m_overwrite == m_overwrite);
      assert(((TripletCalculator*) vc)->m_species[0] == m_species[0]);
      assert(((TripletCalculator*) vc)->m_species[1] == m_species[1]);
      assert(((TripletCalculator*) vc)->m_species[2] == m_species[2]);
      
      MSG_DEBUG("TripletCalculator::setup for " << m_properties.className(), ": registering copy.");    
      
      M_PHASE->registerBondedCalc(vc);
      M_PHASE->registerBondedCalc_0(this);
    }
}

void TripletCalculator::setup()
{
  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

  if(m_symbolName == "undefined")
    throw gError("TripletCalculator::setup", "Attribute 'symbol' undefined!");

  if(m_listName == "undefined")
    throw gError("TripletCalculator::setup", "Attribute 'listName' undefined!");

  
  if(m_species[0] == "undefined")
    throw gError("TripletCalculator::setup", ": Attribute 'species1' has value \"undefined\"!"); 
  if(m_species[1] == "undefined")
    throw gError("TripletCalculator::setup", ": Attribute 'species2' has value \"undefined\"!"); 
  if(m_species[2] == "undefined")
    throw gError("TripletCalculator::setup", ": Attribute 'species3' has value \"undefined\"!"); 
  
  m_firstColour = M_MANAGER->getColour(m_species[0]);
  m_secondColour = M_MANAGER->getColour(m_species[1]);
  m_thirdColour = M_MANAGER->getColour(m_species[2]);
  
  setupTags();

  registerCalc();  
}

void TripletCalculator::setupAfterParticleCreation()
{
  m_boxSize = M_PHASE->boundary()->boundingBox().size();

  // FOLLOWING check cannot be done here because triplet lists 
  // are also created in a setupAfterParticleCreation()

//   // check for consistency of calculator's and triplets' species

//   tripletList* trl = M_PHASE->returnTripletList(m_listName);

//   tripletList::iterator firstTr = trl->begin();

//   if(M_MANAGER->getColour(m_species[0]) != firstTr->a->c)
//     throw gError("TripletCalculator::setupAfterParticleCreation", "Inconsistency between the attribute 'species1' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[1]) != firstTr->b->c)
//     throw gError("TripletCalculator::setupAfterParticleCreation", "Inconsistency between the attribute 'species2' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[2]) != firstTr->c->c)
//     throw gError("TripletCalculator::setupAfterParticleCreation", "Inconsistency between the attribute'species3' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");

}
/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool TripletCalculator::findStage()
{
  
  // currently (2009/07/28) this always returns true
  return Symbol::findStage();
}


// check for consistency of calculator's and triplets' species
void TripletCalculator::checkConsistency() {

  tripletList* trl = M_PHASE->returnTripletList(m_listName);
  
  tripletList::iterator firstTr = trl->begin();
  
  if(M_MANAGER->getColour(m_species[0]) != firstTr->a->c)
    throw gError("TripletCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species1=\"" + m_species[0] + "\" of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  if(M_MANAGER->getColour(m_species[1]) != firstTr->b->c)
    throw gError("TripletCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species2=\"" + m_species[1] + "\"  of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  if(M_MANAGER->getColour(m_species[2]) != firstTr->c->c)
    throw gError("TripletCalculator::findStage", "For module " + name() + ": Inconsistency between the attribute species3=\"" + m_species[2] + "\"  of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");

}   


/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool TripletCalculator::findStage_0()
{
//   // Following check done here because only suitable function after setupAfterParticleCreation() (see also comment there)
//   // check for consistency of calculator's and triplets' species
  
//   tripletList* trl = M_PHASE->returnTripletList(m_listName);
  
//   tripletList::iterator firstTr = trl->begin();
  
//   if(M_MANAGER->getColour(m_species[0]) != firstTr->a->c)
//     throw gError("TripletCalculator::findStage_0", "Inconsistency between the attribute 'species1' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[1]) != firstTr->b->c)
//     throw gError("TripletCalculator::findStage_0", "Inconsistency between the attribute 'species2' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
//   if(M_MANAGER->getColour(m_species[2]) != firstTr->c->c)
//     throw gError("TripletCalculator::findStage_0", "Inconsistency between the attribute'species3' of this calculator and the first species in the assigned list '" + m_listName + "'. This is currently not allowed.");
  
  // currently (2009/07/28) this always returns true
  return Symbol::findStage_0();
}
