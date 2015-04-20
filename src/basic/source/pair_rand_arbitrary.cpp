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



#include "pair_rand_arbitrary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


PairRandArbitrary::PairRandArbitrary(string symbol) :
	ValCalculatorPair(symbol) {
}

PairRandArbitrary::PairRandArbitrary(Simulation* parent) :
	ValCalculatorPair(parent) {
	init();
}

PairRandArbitrary::~PairRandArbitrary() {
}

void PairRandArbitrary::init() {
  m_properties.setClassName("PairRandArbitrary");
  
  STRINGPC
    (symbol, m_symbolName,
     "Name of the random symbol.")
    ;
  
  STRINGPC(pairFactor, m_expression_string,
	   "The mathematical expression to be computed for the pairFactor. It may contain scalars, vectors and matrices, but in the end it must ve of the same type as the final symbol. In the case of vectors or matrices, component-wise multiplication is performed.")
    ;

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Specifies the maximum distance of a pair of particles. Beyond this distance, this module leaves the to be computed symbol zero.")
    ;
  
  STRINGPC
    (useOldFor, m_oldSymbols,
     "Here, you can list the used symbols in the pair-factor, which should be treated as \"old\", i.e., this calculator will not wait for those symbols to be computed beforehand, but it will take what it finds. Separate the symbols by the \"|\"- (\"pipe\"-) symbol."
     )
    ;
  
  INTPC
    (seed, m_seed, 0,
     "Seed to be used for the random number generator."
     )
    ;
  
  m_oldSymbols = "---";
  
  m_cutoff = 0;  

  m_seed = RNG_DEFAULT_SEED;
  
}

//---- Methods ----

void PairRandArbitrary::setup() {
  
  ValCalculatorPair::setup();

  m_rng.setSeed(m_seed);
  
  if (m_cutoff <= 0)
    throw gError("PairRandArbitrary::setup", className() + ": Attribute 'cutoff' has no value > 0!");

  if (m_expression_string == "undefined")
    throw gError("PairRandArbitrary::setup", className() + ": Attribute 'pairFactor' has value \"undefined\"");
  
  MSG_DEBUG("PairRandArbitrary::setup",className() << "  Cutoff used  = " << m_cutoff);

  if (m_allPairs) {
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER,
       cp->setCutoff(m_cutoff);
       cp->setNeedPairs(true);
       vector<ColourPair*>::iterator cpTester = __cp;
       if(++cpTester == __end)
	 {
	   m_function.setExpression(m_expression_string);
	   m_function.setColourPair(cp);
	 }
       )
      ;
  } else /*m_allPairs == false*/
    {
      if (m_species.first == "undefined")
	throw gError("PairRandArbitrary::setup", "For module " + className() + ": Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");
      if (m_species.second == "undefined")
	throw gError("PairRandArbitrary::setup", "For module " + className() + ": Attribute 'species2' has value \"undefined\" and 'allPairs' is disabled.");
      ColourPair* cp= M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);
      cp->setCutoff(m_cutoff);
      cp->setNeedPairs(true);
      m_function.setExpression(m_expression_string);
      m_function.setColourPair(cp);
      
    }
  
  
}

void PairRandArbitrary::setSlot(ColourPair* cp, size_t& slot, bool oneProp) {
  m_slot = slot = cp->tagFormat().addAttribute
    (m_symbolName, m_datatype, false, m_symbolName).offset;
}

bool PairRandArbitrary::findStage() {
  MSG_DEBUG("PairRandArbitrary::findStage", className() << " START: stage = " << m_stage);
  if (m_stage == -1) {
    // no symbols used?
    if (m_function.usedSymbols().empty()) {
      m_stage = 0;
      MSG_DEBUG("PairRandArbitrary::findStage", className() << " for symbol '" << m_symbolName << "': no symbols used: stage is now " << m_stage);
      return true;
    }
    // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
    // this is for setting the stage to '0' when there is no Calculator at all
    // (but probably Integrators or s.th. like that)
    bool nothing = true;
    
    // loop over the used symbols
    
    typed_value_list_t usedSymbols;
    for (typed_value_list_t::const_iterator s = m_function.usedSymbols().begin(); s != m_function.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    
    // go through the string of "old" symbols and remove those from usedSymbols
    if (m_oldSymbols != "---") {
      bool run = true;
      string working = m_oldSymbols;
      MSG_DEBUG("PairRandArbitrary::findStage",className()<<"m_oldSymbols is"<<working);
      while (run) {
	string cur;
	size_t pos = working.find('|');
	
	if (pos == string::npos) {
	  run = false;
	  cur = working;
	} else {
	  cur = string(working, 0, pos);
	  working = string(working, pos+1);
	}
	
	typed_value_list_t toRemove;
	// determine what to remove
	for (typed_value_list_t::const_iterator s = usedSymbols.begin(); s
	       != usedSymbols.end(); ++s) {
	  if ((*s)->name() == cur)
	    toRemove.push_back(*s);
	}
	// remove it
	for (typed_value_list_t::const_iterator s = toRemove.begin(); s
	       != toRemove.end(); ++s)
	  usedSymbols.remove(*s);
	if (toRemove.empty())
	  throw gError("PairRandArbitrary::findStage", "Unable to find old symbol '" + cur + "' among used symbols");
      }
    }
    
    for (typed_value_list_t::const_iterator s = usedSymbols.begin(); s
	   != usedSymbols.end() && !tooEarly; ++s) {
      
      string name = (*s)->name();
      MSG_DEBUG("PairRandArbitrary::findStage", className() << ": now checking for used symbol with complete name " << name);
      // unfortunately, we have to remove the brackets for vectors and tensors
      if (name[0] == '{' || name[0] == '[') {
	// remove the first bracket
	name.erase(0, 1);
	// remove the last bracket; don't know why, but with these arguments it works
	name.erase(name.size()-1, name.size()-1);
      }
      // and the "i", "j" and "ij" of the pair expression have to be removed
      if (name[name.size()-2] == 'i' && name[name.size()-1] == 'j')
	// remove "ij"
	name.erase(name.size()-2, name.size()-1);
      else
	// remove "i" or "j"
	name.erase(name.size()-1, name.size()-1);
      MSG_DEBUG("PairRandArbitrary::findStage", className() << ":short name: " << name);
      
      // we do if((*vCIt) != this) now, so, next should be unnecessary
      /*      if(m_symbolName == name)
	      throw gError("ValCalculatorArbitrary::findStage", className() + " reporting: I cannot use my own symbol as argument.");*/
      
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the
      // stage of it and try to set the stage of this ValCalculator
      // consistently to it,
      
      MSG_DEBUG("PairRandArbitrary::findStage", className() << " for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);
      assert(m_phaseUser == 1 || m_phaseUser == 2);
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
	(
	 M_MANAGER,
	 vector<ValCalculator*>* vCs;
	 /*if(m_phaseUser == 1)*/vCs = &(cp->valCalculatorsFlat());
	 //       if(m_phaseUser == 0) vCs = &(cp->valCalculatorsFlat_0());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin();
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStage", "symbols of VC: ");
		 
		 //             for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 //               cout << *symIt << endl;
		 
		 
		 for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		   {
		     if(*symIt == name)
		       {
			 nothing = false;
			 int stage = (*vCIt)->stage();
			 if(stage == -1)
			   {
			     MSG_DEBUG("PairRandArbitrary::findStage", className() << " for symbol '" << m_symbolName << "': too early because of " << (*vCIt)->className());
			     tooEarly = true;
			     m_stage = -1;
			   }
			 else
			   {
			     if(stage >= m_stage) m_stage = stage+1;
			   }
		       }
		   }
	       }
	   }
	 // can we stop FOR_EACH_COLOUR_PAIR now?
	 if(tooEarly)
	   {
	     __cp = __end;
	     // important because there still comes the ++__cp from the loop
	     --__cp;
	   }
	 )
	;
      
      if (!tooEarly) {
	// we have to search now in the ParticleCaches
	vector<ParticleCache*>* pCs;
	/*if(m_phaseUser == 1)*/pCs
	  = &(Particle::s_cached_flat_properties);
	FOR_EACH
	  (
	   vector<ParticleCache*>,
	   (*pCs)/*Particle::s_cached_flat_properties*/,
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if(*symIt == name)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1)
		     {
		       MSG_DEBUG("PairRandArbitrary::findStage", className() << " for symbol '" << m_symbolName << "': too early because of " << (*__iFE)->className());
		       tooEarly = true;
		       m_stage = -1;
		     }
		   else
		     {
		       if(stage >= m_stage) m_stage = stage+1;
		     }
		 }
	     }
	   // may we abort the loop over the ParticleCaches?
	   if(tooEarly)
	     {
	       __iFE = __end;
	       // important because there still comes the ++__iFE from the loop
	       --__iFE;
	     }
	   );
	
      }
    } // end loop over symbols
    if (tooEarly)
      return false;
    if (m_stage == -1) {
      if (nothing) {
	m_stage = 0;
	MSG_DEBUG("PairRandArbitrary::findStage", className() << " for symbol '" << m_symbolName << "': stage is now " << m_stage);
	return true;
      } else
	return false;
    } else {
      MSG_DEBUG("PairRandArbitrary::findStage", className() << " for symbol '" << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else {
    MSG_DEBUG("PairRandArbitrary::findStage", className() << " for symbol '" << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}

bool PairRandArbitrary::findStage_0() {
  MSG_DEBUG("PairRandArbitrary::findStage_0", className() << " START: stage = " << m_stage);
  if (m_stage == -1) {
    // no symbols used?
    if (m_function.usedSymbols().empty()) {
      m_stage = 0;
      MSG_DEBUG("PairRandArbitrary::findStage_0", className() << " for symbol '" << m_symbolName << "': no symbols used: stage is now " << m_stage);
      return true;
    }
    // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
    // this is for setting the stage to '0' when there is no Calculator at all
    // (but probably Integrators or s.th. like that)
    bool nothing = true;
    
    // loop over the used symbols
    
    typed_value_list_t usedSymbols;
    for (typed_value_list_t::const_iterator s = m_function.usedSymbols().begin(); s != m_function.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    // go through the string of "old" symbols and remove those from usedSymbols
    if (m_oldSymbols != "---") {
      bool run = true;
      string working = m_oldSymbols;
      while (run) {
	string cur;
	size_t pos = working.find('|');
	
	if (pos == string::npos) {
	  run = false;
	  cur = working;
	} else {
	  cur = string(working, 0, pos);
	  working = string(working, pos+1);
	}
	
	typed_value_list_t toRemove;
	// determine what to remove
	for (typed_value_list_t::const_iterator s = usedSymbols.begin(); s
	       != usedSymbols.end(); ++s) {
	  if ((*s)->name() == cur)
	    toRemove.push_back(*s);
	}
	// remove it
	for (typed_value_list_t::const_iterator s = toRemove.begin(); s
	       != toRemove.end(); ++s)
	  usedSymbols.remove(*s);
	if (toRemove.empty())
	  throw gError("PairRandArbitrary::findStage_0", "Unable to find old symbol '" + cur + "' among used symbols");
      }
    }
    
    for (typed_value_list_t::const_iterator s = usedSymbols.begin(); s
	   != usedSymbols.end() && !tooEarly; ++s) {
      
      string name = (*s)->name();
      MSG_DEBUG("PairRandArbitrary::findStage_0", className() << ": now checking for used symbol with complete name " << name);
      // unfortunately, we have to remove the brackets for vectors and tensors
      if (name[0] == '{' || name[0] == '[') {
	// remove the first bracket
	name.erase(0, 1);
	// remove the last bracket; don't know why, but with these arguments it works
	name.erase(name.size()-1, name.size()-1);
      }
      // and the "i", "j" and "ij" of the pair expression have to be removed
      if (name[name.size()-2] == 'i' && name[name.size()-1] == 'j')
	// remove "ij"
	name.erase(name.size()-2, name.size()-1);
      else
	// remove "i" or "j"
	name.erase(name.size()-1, name.size()-1);
      MSG_DEBUG("PairRandArbitrary::findStage_0", className() << ":short name: " << name);
      
      // we do if((*vCIt) != this) now, so, next should be unnecessary
      /*      if(m_symbolName == name)
	      throw gError("ValCalculatorArbitrary::findStage_0", className() + " reporting: I cannot use my own symbol as argument.");*/
      
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the
      // stage of it and try to set the stage of this ValCalculator
      // consistently to it,
      
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
	(
	 M_MANAGER,
	 assert(m_phaseUser == 2 || m_phaseUser == 0);
	 vector<ValCalculator*>* vCs;
	 //       if(m_phaseUser == 1) vCs = &(cp->valCalculatorsFlat());
	 /*if(m_phaseUser == 0)*/vCs = &(cp->valCalculatorsFlat_0());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin();
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStage_0", "symbols of VC: ");
		 
		 //             for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 //               cout << *symIt << endl;
		 
		 
		 for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		   {
		     if(*symIt == name)
		       {
			 nothing = false;
			 int stage = (*vCIt)->stage();
			 if(stage == -1)
			   {
			     MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '" << m_symbolName << "': too early because of " << (*vCIt)->className());
			     tooEarly = true;
			     m_stage = -1;
			   }
			 else
			   {
			     if(stage >= m_stage) m_stage = stage+1;
			   }
		       }
		   }
	       }
	   }
	 // can we stop FOR_EACH_COLOUR_PAIR now?
	 if(tooEarly)
	   {
	     __cp = __end;
	     // important because there still comes the ++__cp from the loop
	     --__cp;
	   }
	 );
      
      if (!tooEarly) {
	// we have to search now in the ParticleCaches
	vector<ParticleCache*>* pCs;
	//             if(m_phaseUser == 1) pCs = &(Particle::s_cached_flat_properties);
	/*if(m_phaseUser == 0)*/pCs
	  = &(Particle::s_cached_flat_properties_0);
	//             MSG_DEBUG("ValCalculatorArbitrary::findStage_0", "loop over ParticleCaches START");
	FOR_EACH
	  (
	   vector<ParticleCache*>,
	   (*pCs)/*Particle::s_cached_flat_properties*/,
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if(*symIt == name)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1)
		     {
		       MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '" << m_symbolName << "': too early because of " << (*__iFE)->className());
		       tooEarly = true;
		       m_stage = -1;
		     }
		   else
		     {
		       if(stage >= m_stage) m_stage = stage+1;
		     }
		 }
	     }
	   // may we abort the loop over the ParticleCaches?
	   if(tooEarly)
	     {
	       __iFE = __end;
	       // important because there still comes the ++__iFE from the loop
	       --__iFE;
	     }
	   );
	
      }
    } // end loop over symbols
    if (tooEarly)
      return false;
    if (m_stage == -1) {
      if (nothing) {
	m_stage = 0;
	MSG_DEBUG("PairRandArbitrary::findStage_0", className() << " for symbol '" << m_symbolName << "': stage is now " << m_stage);
	return true;
      } else
	return false;
    } else {
      MSG_DEBUG("PairRandArbitrary::findStage_0", className() << " for symbol '" << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else {
    MSG_DEBUG("PairRandArbitrary::findStage_0", className() << " for symbol '" << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}


