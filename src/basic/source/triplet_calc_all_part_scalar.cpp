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
#include "triplet_calc_all_part_scalar.h"


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

const SymbolRegister<TripletCalcAllPartScalar> triplet_calc_all_part_scalar("TripletCalcAllPartScalar");

TripletCalcAllPartScalar::TripletCalcAllPartScalar(Simulation *simulation): TripletCalcPartScalar(simulation)
{  
  init();
}


void TripletCalcAllPartScalar::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalcAllPartScalar");

  m_properties.setDescription("TripletCalculator for cached properties computed during a loop over bonded triplets. This \a Symbol computes (possibly different) runtme-compiled expressions depending on the cosine \"ca\" of the triplet angle and sums them up in the respective particles of the triplet. By default, attribute 'expression' denotes the value computed and stored for each particle in the triplet. With attributes 'expression1st' and 'expression3rd' this can be modified for the two outer (1st, 3rd) particles. By default, attribute 'symbol' denotes the variable for storing the value for each particle in the triplet. With attributes 'symbol1st' and 'symbol3rd' this can be modified for the two outer (1st, 3rd) particles. ");
  
  FUNCTIONFIXEDPC
      (expression1st, m_expression1st, "Scalar function depending on the variable \"ca\" which is the cosine of the triplet-angle alpha. If defined, it denotes the value computed and stored in the first particle of the triplet. Otherwise the definition of attribute 'expression' is used.");

  m_expression1st.addVariable("ca");
  m_expression1st.setExpression("Same as attribute 'expression'.");

  FUNCTIONFIXEDPC
      (expression3rd, m_expression3rd, "Scalar function depending on the variable \"ca\" which is the cosine of the triplet-angle alpha. If defined, it denotes the value computed and stored in the 3rd particle of the triplet. Otherwise the definition of attribute 'expression' is used.");

  m_expression3rd.addVariable("ca");
  m_expression3rd.setExpression("Same as attribute 'expression'.");

  STRINGPC
      (symbol1st, m_symbolName1st, "Symbol name for storing the result in the 1st particle of the triplet. If not defined, the name defined in attribute 'symbol' is taken.");

  m_symbolName1st = "Same as attribute 'symbol'";

  STRINGPC
      (symbol3rd, m_symbolName3rd, "Symbol name for storing the result in the 3rd particle of the triplet. If not defined, the name defined in attribute 'symbol' is taken.");

  m_symbolName3rd = "Same as attribute 'symbol'";

}


void TripletCalcAllPartScalar::setup()
{

  // works because we redefine setupTags()
  TripletCalcPartScalar::setup();

  if(m_expression1st.expression() == "Same as attribute 'expression'.")
    m_expression1st.setExpression(m_expression.expression());

  if(m_expression3rd.expression() == "Same as attribute 'expression'.")
    m_expression3rd.setExpression(m_expression.expression());

  if(m_symbolName1st == "Same as attribute 'symbol'.")
    m_symbolName1st = m_symbolName;

  if(m_symbolName3rd == "Same as attribute 'symbol'.")
    m_symbolName3rd = m_symbolName;

}


void TripletCalcAllPartScalar::setupTags()
{

  if(m_overwrite)
    {
      
      try
	{
	  m_slots[0] = Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].indexOf(m_symbolName1st, m_datatype);
	  m_slots[0] = Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].offsetByIndex(m_slots[0]);
	}
      catch(gError& err)
	{
	  throw gError("TripletCalcAllPartScalar::setup", ": search for symbol for species '" + m_species[1] + " failed. The message was " + err.message()); 
	}
      
      try
	{
	  m_slots[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].indexOf(m_symbolName, m_datatype);
	  m_slots[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].offsetByIndex(m_slots[1]);
	}
      catch(gError& err)
	{
	  throw gError("TripletCalcAllPartScalar::setup", ": search for symbol for species '" + m_species[1] + " failed. The message was " + err.message()); 
	}
      
      try
	{
	  m_slots[2] = Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].indexOf(m_symbolName3rd, m_datatype);
	  m_slots[2] = Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].offsetByIndex(m_slots[2]);
	}
      catch(gError& err)
	{
	  throw gError("TripletCalcAllPartScalar::setup", ": search for symbol for species '" + m_species[2] + " failed. The message was " + err.message()); 
	}

    }
  else // m_overwrite = false
    {
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].attrExists(m_symbolName1st)) 
        throw gError("TripletCalcAllPartScalar::setup", ": Symbol " + m_symbolName1st + " is already existing for species '" + m_species[0] + "'. Second definition is not allowed for overwrite = \"no\"");
      
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].attrExists(m_symbolName)) 
        throw gError("TripletCalcAllPartScalar::setup", ": Symbol " + m_symbolName + " is already existing for species '" + m_species[1] + "'. Second definition is not allowed for overwrite = \"no\"");
      
      if(Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].attrExists(m_symbolName3rd)) 
        throw gError("TripletCalcAllPartScalar::setup", ": Symbol " + m_symbolName3rd + " is already existing for species '" + m_species[2] + "'. Second definition is not allowed for overwrite = \"no\"");
      
      // see CONVENTION5 for rule about persistencies
      m_slots[0] = Particle::s_tag_format[M_MANAGER->getColour(m_species[0])].addAttribute(m_symbolName1st, m_datatype, /*persist.first*/false, m_symbolName).offset;
      
      if(m_species[0] != m_species[1] || m_symbolName1st != m_symbolName)
	{
	  // see CONVENTION5 for rule about persistencies
	  m_slots[1] = Particle::s_tag_format[M_MANAGER->getColour(m_species[1])].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
	}
      else m_slots[1] = m_slots[0];
      
      if(m_species[2] == m_species[0] && m_symbolName3rd == m_symbolName1st)
	m_slots[2] = m_slots[0];
      else {
	if(m_species[2] == m_species[1] && m_symbolName3rd == m_symbolName)
	  m_slots[2] = m_slots[1];
	else
	  // see CONVENTION5 for rule about persistencies
	  m_slots[2] = Particle::s_tag_format[M_MANAGER->getColour(m_species[2])].addAttribute(m_symbolName3rd, m_datatype, /*persist.first*/false, m_symbolName).offset;
      }
      
    } // end of else of if(m_overwrite == false)        

}


void TripletCalcAllPartScalar::setupAfterParticleCreation()
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
bool TripletCalcAllPartScalar::findStage()
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
bool TripletCalcAllPartScalar::findStage_0()
{
  // currently (2010/05/17) this always returns true
  return TripletCalcPartScalar::findStage_0();
}

#ifdef _OPENMP
// FIXME: This module is not yet parallelised. If you parallelise, the following function could roughly do what is commented out now
void TripletCalcAllPartScalar::mergeCopies(size_t thread_no) {
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
