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
#include "triplet_calc_central_part_scalar.h"

// #include "triplet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

const SymbolRegister<TripletCalcCentralPartScalar> triplet_calc_central_part_scalar("TripletCalcCentralPartScalar");

TripletCalcCentralPartScalar::TripletCalcCentralPartScalar(Simulation *simulation): TripletCalcPartScalar(simulation)
{  
  init();
}


void TripletCalcCentralPartScalar::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalcCentralPartScalar");

  m_properties.setDescription("TripletCalculator for cached properties computed during a loop over bonded triplets. This TripletCalculator specifically computes a scalar user-defined property depending on the cosine of the triplet angle \"ca\". For each triplet this property only contributes to the CENTRE particle by summation! ");
  
}


void TripletCalcCentralPartScalar::setup()
{
  // Since registration of Symbols is only necessary for the central particle, 
  // we don't do the following line and instead perform all the setup here.
  // The reason is that TripletCalcPartScalar::setup calls TripletCalculator::setup,
  // which sets all the symbols and their slots as if they had the same physical meaning, 
  // which some inheriting classes either do not require (as here) or find to restrictive.

  //  TripletCalcPartScalar::setup();

  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

  if(m_symbolName == "undefined")
    throw gError("TripletCalcCentralPartScalar::setup", "Attribute 'symbol' undefined!");

  if(m_listName == "undefined")
    throw gError("TripletCalcCentralPartScalar::setup", "Attribute 'listName' undefined!");

    if(m_species[0] == "undefined")
      throw gError("TripletCalcCentralPartScalar::setup", ": Attribute 'species1' has value \"undefined\"!"); 
    if(m_species[1] == "undefined")
      throw gError("TripletCalcCentralPartScalar::setup", ": Attribute 'species2' has value \"undefined\"!"); 
    if(m_species[2] == "undefined")
      throw gError("TripletCalcCentralPartScalar::setup", ": Attribute 'species3' has value \"undefined\"!"); 

    m_firstColour = M_MANAGER->getColour(m_species[0]);
    m_secondColour = M_MANAGER->getColour(m_species[1]);
    m_thirdColour = M_MANAGER->getColour(m_species[2]);


    // add symbols only to species of central particle

    if(m_overwrite)
    {
      try {
	// first the index
	m_slots[1] = Particle::s_tag_format[m_secondColour].indexOf(m_symbolName, m_datatype);
	// now the memory-offset from the index
	m_slots[1] = Particle::s_tag_format[m_secondColour].offsetByIndex(m_slots[1]);
      }
      catch(gError& err) {
	throw gError("TripletCalcCentralPartScalar::setup", ": search for symbol for species '" + m_species[1] + " failed. The message was " + err.message()); 
      }
      
    }
    else // m_overwrite = false
    {
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].attrExists(m_symbolName)) 
        throw gError("TripletCalcCentralPartScalar::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[1] + "'. Second definition is not allowed for overwrite = \"no\"");

      MSG_DEBUG("TripletCalcCentralPartScalar::setup", "adding symbol with name " << m_symbolName);
      
      // see CONVENTION5 for rule about persistencies
      m_slots[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;

    } // end of else of if(m_overwrite == false)        
    
    if(m_phaseUser == 0)    
      M_PHASE->registerBondedCalc_0(this);
    else if(m_phaseUser == 1)    
      M_PHASE->registerBondedCalc(this);
    else // so it is 2
    {
      TripletCalculator* vc = copyMySelf();
      assert(((TripletCalcCentralPartScalar*) vc)->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((TripletCalcCentralPartScalar*) vc)->m_slots[0] == m_slots[0]);
      assert(((TripletCalcCentralPartScalar*) vc)->m_slots[1] == m_slots[1]);
      assert(((TripletCalcCentralPartScalar*) vc)->m_slots[2] == m_slots[2]);
      assert(((TripletCalcCentralPartScalar*) vc)->m_parent == m_parent);
      assert(((TripletCalcCentralPartScalar*) vc)->m_overwrite == m_overwrite);
      assert(((TripletCalcCentralPartScalar*) vc)->m_species[0] == m_species[0]);
      assert(((TripletCalcCentralPartScalar*) vc)->m_species[1] == m_species[1]);
      assert(((TripletCalcCentralPartScalar*) vc)->m_species[2] == m_species[2]);
        
      MSG_DEBUG("TripletCalcCentralPartScalar::setup", ": registering copy.");    
       
      M_PHASE->registerBondedCalc(vc);
      M_PHASE->registerBondedCalc_0(this);
    }

}

void TripletCalcCentralPartScalar::setupAfterParticleCreation()
{
  TripletCalcPartScalar::setupAfterParticleCreation();
}

/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool TripletCalcCentralPartScalar::findStage()
{
  
  // currently (2010/05/17) this always returns true
  return TripletCalcPartScalar::findStage();
}

/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool TripletCalcCentralPartScalar::findStage_0()
{
  // currently (2010/05/17) this always returns true
  return TripletCalcPartScalar::findStage_0();
}

#ifdef _OPENMP
// FIXME: This module is not yet parallelised. If you parallelise, the following function could roughly do what is commented out now
void TripletCalcCentralPartScalar::mergeCopies(size_t thread_no) {
//   size_t slot2 = m_slots[1];

//   size_t copySlot2 = m_copy_slots[thread_no][1];

//   size_t vecSlot2 = m_vector_slots[1];

//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, m_secondColour,
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.doubleByOffset(slot2)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j] = 0;
//     }
//   );
}
#endif
