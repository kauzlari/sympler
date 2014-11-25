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


#include "bonded_pair_particle_calc.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


BondedPairParticleCalc::BondedPairParticleCalc(Simulation* parent)
  : ValCalculatorPart(parent)
{
  init();
}

BondedPairParticleCalc::~BondedPairParticleCalc()
{
}

void BondedPairParticleCalc::init()
{
//   m_properties.setClassName("BondedPairParticleCalc");
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("BondedPairParticleCalc");

  STRINGPC(listName, m_listName, "Identifier of the list of bonded pairs, this Calculator should belong to.");

  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the calculated property.");


  BOOLPC
      (overwrite, m_overwrite,
       "Is this calculator allowed to overwrite already existing symbols " 
           "with name 'symbol' ?");

  INTPC
      (symmetry, m_symmetry, -2,
       "How does the pair expression behave under index interchange? \"1\" means symmetric behaviour. \"-1\" means antisymmetric behaviour.");


#ifdef _OPENMP
  m_particleCalculator = true;
#endif
  m_overwrite = false;
  m_listName = "undefined";

  m_symmetry = -2;

}

void BondedPairParticleCalc::setup()
{
   if(m_symmetry != -1 && m_symmetry != 1)
     throw gError("BondedPairParticleCalc::setup", "For module " + name() + ": Define attribute 'symmetry'! It may only be \"-1\" (antisymmetry) or \"1\" (symmetry). Currently it is " + ObjToString(m_symmetry) + "!");

  if(m_listName == "undefined")
    throw gError("BondedPairParticleCalc::setup", "Attribute 'listName' is undefined!");
       
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("BondedPairParticleCalc::setup", ": Attribute 'stage' is " + ObjToString(m_phaseUser) + ", which is none of the allowed values \"0\", \"1\", \"2\".");
   
  if(m_species.first == "undefined")
    throw gError("BondedPairParticleCalc::setup", "For module " + name() + ": Attribute 'species1' has value \"undefined\"."); 
  if(m_species.second == "undefined")
    throw gError("BondedPairParticleCalc::setup", "For module " + name() + ": Attribute 'species2' has value \"undefined\"."); 
  
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);
  
  
  // next lines are just necessary for non-bonded pairs, so not here!
  //      cp->setNeedPairs(true);
  //   cp->setCutoff(m_cutoff);
  
  if(m_overwrite)
    {
      try
	{
	  m_slots.first = Particle::s_tag_format[cp->firstColour()].indexOf(m_symbolName, m_datatype);
	  m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByIndex(m_slots.first);
	}
      catch(gError& err)
	{
	  throw gError("BondedPairParticleCalc::setup", ": search for symbol for species '" + M_MANAGER->species(cp->firstColour()) + " failed, but need one because 'overwrite = \"yes\"'. The message was " + err.message()); 
	}
      try
	{
	  m_slots.second = Particle::s_tag_format[cp->secondColour()].indexOf(m_symbolName, m_datatype);
	  m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByIndex(m_slots.second);
	}
      catch(gError& err)
	{
	  throw gError("BondedPairParticleCalc::setup", ": search for symbol for species '" + M_MANAGER->species(cp->secondColour()) + " failed, but need one because 'overwrite = \"yes\"'. The message was " + err.message()); 
	}      
    }
  else // m_overwrite = false
    {
      if(Particle::s_tag_format[cp->firstColour()].attrExists(m_symbolName)) 
        throw gError("BondedPairParticleCalc::setup", ": Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(cp->firstColour()) + "'. Second definition is not allowed for overwrite = \"no\"");
      
      if(Particle::s_tag_format[cp->secondColour()].attrExists(m_symbolName))
        throw gError("BondedPairParticleCalc::setup", ": Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(cp->secondColour()) + "'. Second definition is not allowed for overwrite = \"no\"");
      
      // see CONVENTION5 for rule about persistencies
      m_slots.first = Particle::s_tag_format[cp->firstColour()].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
      
      if(cp->firstColour() != cp->secondColour())
	{
	  // see CONVENTION5 for rule about persistencies
	  m_slots.second = Particle::s_tag_format[cp->secondColour()].addAttribute(m_symbolName, m_datatype, /*persist.second*/false, m_symbolName).offset;
	}
      else m_slots.second = m_slots.first;
    } // end of else-case where m_overwrite = false
  
  MSG_DEBUG("BondedPairParticleCalc::setup", ": registering " << this->name() << ".");
  
  if(m_phaseUser == 0)    
    cp->registerBondedCalc_0(this);
  else if(m_phaseUser == 1)    
    cp->registerBondedCalc(this);
  else // so it is 2
    {
      ValCalculator* vc = copyMySelf() /*new BondedPairParticleCalc(*this)*/;

      copyMembersTo(vc);
      
      MSG_DEBUG("BondedPairParticleCalc::setup", ": registering copy of " << this->name() << ", CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");    
      
      cp->registerBondedCalc(vc);
      cp->registerBondedCalc_0(this);
    }

}

void BondedPairParticleCalc::copyMembersTo(ValCalculator* vc)
{

  ((BondedPairParticleCalc*) vc)->m_overwrite = m_overwrite;
  ((BondedPairParticleCalc*) vc)->m_symmetry = m_symmetry;
#ifdef _OPENMP
  ((BondedPairParticleCalc*) vc)->m_particleCalculator = m_particleCalculator;
#endif

}


#ifdef _OPENMP
void BondedPairParticleCalc::mergeCopies(ColourPair* cp, int thread_no) {

  // this function currently does nothing because the parallel version computes as the serial version and writes directly in the main particle slots. The following code is from pair_particle_vector.cpp and shows how it could work

//   size_t slot1 = m_slots.first;
//   size_t slot2 = m_slots.second;

//   size_t copySlot1 = m_copy_slots[thread_no].first;
//   size_t copySlot2 = m_copy_slots[thread_no].second;
//   size_t vecSlot1 = m_vector_slots.first;
//   size_t vecSlot2 = m_vector_slots.second;

//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, cp->firstColour(),
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       i->tag.pointByOffset(slot1)[j] += (*i->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j];
//       (*i->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j] = 0;
//     }
//   );
//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, cp->secondColour(),
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       i->tag.pointByOffset(slot2)[j] += (*i->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j];
//       (*i->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j] = 0;
//     }
//   );
}
#endif


