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



#include "val_calculator_rho.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

const SymbolRegister<ValCalculatorRho> val_calc_rho("ValCalculatorRho");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

ValCalculatorRho::ValCalculatorRho(string symbol, WeightingFunction *wf)
  : NonBondedPairParticleCalculator(symbol), m_cutoff(wf->cutoff()), m_wf(wf)
{}

ValCalculatorRho::ValCalculatorRho(/*Node*/Simulation* parent)
  : NonBondedPairParticleCalculator(parent)
{
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  init();
}

void ValCalculatorRho::init()
{
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("ValCalculatorRho");

  m_properties.setDescription("Calculates the local density via pair summation.");

  STRINGPC
      (symbol, m_symbolName,
        "Name of the symbol for the calculated density.");

  STRINGPC
      (weightingFunction, m_wfName,
        "Name of the weighting function to be used.");

  BOOLPC
      (selfContribution, m_self, "Should the self contribution of the particles be included into the local density calculation?");

  m_self = true;
  m_wfName = "default";

#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}

void ValCalculatorRho::setup()
{
  m_wf = M_SIMULATION->findWeightingFunction(m_wfName);

  m_cutoff = m_wf -> cutoff();

  if(m_allPairs) {
    
    size_t colour;
    // is the symbol already existing somewhere?
    for (colour = 0; colour < M_MANAGER->nColours(); ++colour) {
      if(Particle::s_tag_format[colour].attrExists(m_symbolName))
        throw gError("ValCalculatorRho::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(colour) + "'. Second definition is not allowed with this Symbol calculator");
    }

    colour = 0;

    // see CONVENTION5 for rule about persistency
    m_slots.first = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persistency*/false, m_symbolName).offset;
    m_slots.second = m_slots.first;
    // ValCalculatorRho should create PCADensitySelf if self contribution is switched on
    if(m_self) {
      if(m_phaseUser == 0)
        Particle::registerCache_0
	  (new ParticleCacheDensitySelfContribution
	   (colour, m_slots.first, m_wf, m_symbolName));
      if(m_phaseUser == 1)
        Particle::registerCache
	  (new ParticleCacheDensitySelfContribution
	   (colour, m_slots.first, m_wf, m_symbolName));
      if(m_phaseUser == 2) {
	Particle::registerCache_0
	  (new ParticleCacheDensitySelfContribution
	   (colour, m_slots.first, m_wf, m_symbolName));
        Particle::registerCache
	  (new ParticleCacheDensitySelfContribution
	   (colour, m_slots.first, m_wf, m_symbolName));
      }
    }
    
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER,
       
       cp->setCutoff(m_cutoff);
       cp->setNeedPairs(true);
       
       if(colour < cp->firstColour() || colour < cp->secondColour()) {
	 ++colour;
	 if(cp->firstColour() == colour) {
	   // see CONVENTION5 for rule about persistencies
	   m_slots.first = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persistency*/false, m_symbolName).offset;
	   // ValCalculatorRho should create PCADensitySelf if self contribution is switched on
	   if(m_self) {
	     if(m_phaseUser == 0)
	       Particle::registerCache_0
		 (new ParticleCacheDensitySelfContribution
                  (colour, m_slots.first, m_wf, m_symbolName));
	     if(m_phaseUser == 1)
	       Particle::registerCache
		 (new ParticleCacheDensitySelfContribution
                  (colour, m_slots.first, m_wf, m_symbolName));
	     if(m_phaseUser == 2) {
	       Particle::registerCache_0
		 (new ParticleCacheDensitySelfContribution
                  (colour, m_slots.first, m_wf, m_symbolName));
	       Particle::registerCache
		 (new ParticleCacheDensitySelfContribution
                  (colour, m_slots.first, m_wf, m_symbolName));
	     }
	   }
	 }
	 else // here we have to make sure at least that the offsets are correct
	   m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByName(m_symbolName);

	 if(cp->secondColour() == colour) {
	   // see CONVENTION5 for rule about persistencies
	   m_slots.second = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persistency*/false, m_symbolName).offset;
	   // ValCalculatorRho should create PCADensitySelf if self contribution is switched on
	   if(m_self) {
	     if(m_phaseUser == 0)
	       Particle::registerCache_0
		 (new ParticleCacheDensitySelfContribution
                  (colour, m_slots.second, m_wf, m_symbolName));
	     if(m_phaseUser == 1)
	       Particle::registerCache
		 (new ParticleCacheDensitySelfContribution
		  (colour, m_slots.second, m_wf, m_symbolName));
	     if(m_phaseUser == 2) {
	       Particle::registerCache_0
		 (new ParticleCacheDensitySelfContribution
		  (colour, m_slots.second, m_wf, m_symbolName));
	       Particle::registerCache
		 (new ParticleCacheDensitySelfContribution
		  (colour, m_slots.second, m_wf, m_symbolName));
	     }
	   }
	 }
	 else // here we have to make sure at least that the offsets are correct
	   m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByName(m_symbolName);

       }  // end of if(colour < cp->firstColour() || colour < cp->secondColour())
       else {
	 // here we have to make sure at least that the offsets are correct
	 m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByName(m_symbolName);
	 m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByName(m_symbolName);
       }
       
       vector<ColourPair*>::iterator cpTester = __cp;
       // is it the last calculator to be created?
       if(++cpTester == __end) {
	 if(m_phaseUser == 0)
	   cp->registerCalc_0(this);
	 else if(m_phaseUser == 1)
	   cp->registerCalc(this);
	 else { // so it is 2 
	   ValCalculator* vc = /*copyMySelf()*/new ValCalculatorRho(*this);
	   assert(((ValCalculatorRho*) vc)->m_symbolName == m_symbolName);
	   assert(vc->stage() == m_stage);
	   assert(((ValCalculatorRho*) vc)->m_wf == m_wf);
	   assert(((ValCalculatorRho*) vc)->m_wfName == m_wfName);
	   assert(((ValCalculatorRho*) vc)->m_slots.first == m_slots.first);
	   assert(((ValCalculatorRho*) vc)->m_slots.second == m_slots.second);
	   assert(((ValCalculatorRho*) vc)->m_parent == m_parent);
	   assert(((ValCalculatorRho*) vc)->m_species.first == m_species.first);
	   assert(((ValCalculatorRho*) vc)->m_species.second == m_species.second);
	   
	   cp->registerCalc(vc);
	   cp->registerCalc_0(this);
	 }
       }
       // No? Then make a copy
       else {
	 ValCalculator* vc = /*copyMySelf()*/new ValCalculatorRho(*this);
	 assert(((ValCalculatorRho*) vc)->m_symbolName == m_symbolName);
	 assert(vc->stage() == m_stage);
	 assert(((ValCalculatorRho*) vc)->m_wf == m_wf);
	 assert(((ValCalculatorRho*) vc)->m_wfName == m_wfName);
	 assert(((ValCalculatorRho*) vc)->m_slots.first == m_slots.first);
	 assert(((ValCalculatorRho*) vc)->m_slots.second == m_slots.second);
	 assert(((ValCalculatorRho*) vc)->m_parent == m_parent);
	 assert(((ValCalculatorRho*) vc)->m_species.first == m_species.first);
	 assert(((ValCalculatorRho*) vc)->m_species.second == m_species.second);
	 
	 if(m_phaseUser == 0)
	   cp->registerCalc_0(vc);
	 else if(m_phaseUser == 1)
	   cp->registerCalc(vc);
	 else { // so it is 2
	   cp->registerCalc(vc);
	   
	   vc = /*copyMySelf()*/new ValCalculatorRho(*this);
	   assert(((ValCalculatorRho*) vc)->m_symbolName == m_symbolName);
	   assert(vc->stage() == m_stage);
	   assert(((ValCalculatorRho*) vc)->m_wf == m_wf);
	   assert(((ValCalculatorRho*) vc)->m_wfName == m_wfName);
	   assert(((ValCalculatorRho*) vc)->m_slots.first == m_slots.first);
	   assert(((ValCalculatorRho*) vc)->m_slots.second == m_slots.second);
	   assert(((ValCalculatorRho*) vc)->m_parent == m_parent);
	   assert(((ValCalculatorRho*) vc)->m_species.first == m_species.first);
	   assert(((ValCalculatorRho*) vc)->m_species.second == m_species.second);
	   
	   cp->registerCalc_0(vc);
	 }
       }
       );
  }
  else { /*m_allPairs == false case*/
    if(m_species.first == "undefined")
      throw gError("ValCalculatorRho::setup", "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");
    if(m_species.second == "undefined")
      throw gError("ValCalculatorRho::setup", "Attribute 'species2' has value \"undefined\" and 'allPairs' is disabled.");
    
    ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
    
    cp->setCutoff(m_cutoff);
    cp->setNeedPairs(true);
    
    if(Particle::s_tag_format[cp->firstColour()].attrExists(m_symbolName))
      throw gError("ValCalculatorRho::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(cp->firstColour()) + "'. Second definition is not allowed with this Symbol calculator");
    if(Particle::s_tag_format[cp->secondColour()].attrExists(m_symbolName))
      throw gError("ValCalculatorRho::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(cp->secondColour()) + "'. Second definition is not allowed with this Symbol calculator");
    
    // see CONVENTION5 for rule about persistencies
    m_slots.first = Particle::s_tag_format[cp->firstColour()].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
    if(m_self) {
      if(m_phaseUser == 0)
        Particle::registerCache_0
	  (new ParticleCacheDensitySelfContribution
	   (cp->firstColour(), m_slots.first, m_wf, m_symbolName));
      if(m_phaseUser == 1)
        Particle::registerCache
	  (new ParticleCacheDensitySelfContribution
	   (cp->firstColour(), m_slots.first, m_wf, m_symbolName));
      if(m_phaseUser == 2) {
        Particle::registerCache_0
	  (new ParticleCacheDensitySelfContribution
	   (cp->firstColour(), m_slots.first, m_wf, m_symbolName));
        Particle::registerCache
	  (new ParticleCacheDensitySelfContribution
	   (cp->firstColour(), m_slots.first, m_wf, m_symbolName));
      }
    }
    
    if(cp->firstColour() != cp->secondColour())
      {
	// see CONVENTION5 for rule about persistencies
	m_slots.second = Particle::s_tag_format[cp->secondColour()].addAttribute(m_symbolName, m_datatype, /*persistency*/false, m_symbolName).offset;
	if(m_self) {
	  if(m_phaseUser == 0)
	    Particle::registerCache_0
              (new ParticleCacheDensitySelfContribution
	       (cp->secondColour(), m_slots.second, m_wf, m_symbolName));
	  if(m_phaseUser == 1)
	    Particle::registerCache
              (new ParticleCacheDensitySelfContribution
	       (cp->secondColour(), m_slots.second, m_wf, m_symbolName));
	  if(m_phaseUser == 2) {
	    Particle::registerCache_0
              (new ParticleCacheDensitySelfContribution
	       (cp->secondColour(), m_slots.second, m_wf, m_symbolName));
	    Particle::registerCache
              (new ParticleCacheDensitySelfContribution
	       (cp->secondColour(), m_slots.second, m_wf, m_symbolName));
	  }
	}
      }
    else m_slots.second = m_slots.first;
    
    if(m_phaseUser == 0)
      cp->registerCalc_0(this);
    if(m_phaseUser == 1)
      cp->registerCalc(this);
    else { // so it is 2
      ValCalculator* vc = /*copyMySelf()*/new ValCalculatorRho(*this);
      assert(((ValCalculatorRho*) vc)->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((ValCalculatorRho*) vc)->m_wf == m_wf);
      assert(((ValCalculatorRho*) vc)->m_wfName == m_wfName);
      assert(((ValCalculatorRho*) vc)->m_slots.first == m_slots.first);
      assert(((ValCalculatorRho*) vc)->m_slots.second == m_slots.second);
      assert(((ValCalculatorRho*) vc)->m_parent == m_parent);
      assert(((ValCalculatorRho*) vc)->m_species.first == m_species.first);
      assert(((ValCalculatorRho*) vc)->m_species.second == m_species.second);
      
      cp->registerCalc_0(this);
      cp->registerCalc(vc);
    }
  }
}


#ifdef _OPENMP
void ValCalculatorRho::mergeCopies(ColourPair* cp, int thread_no) {
  size_t slot1 = m_slots.first;
  size_t slot2 = m_slots.second;

  size_t copySlot1 = m_copy_slots[thread_no].first;
  size_t copySlot2 = m_copy_slots[thread_no].second;
  size_t vecSlot1 = m_vector_slots.first;
  size_t vecSlot2 = m_vector_slots.second;

  FOR_EACH_PARTICLE_C
  (M_PHASE, cp->firstColour(),
      __iSLFE->tag.doubleByOffset(slot1) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1];
        (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1] = 0;
  );

  // second colour!
  FOR_EACH_PARTICLE_C
  (M_PHASE, cp->secondColour(),
      __iSLFE->tag.doubleByOffset(slot2) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2];
        (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2] = 0;
  );
}
#endif
