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


#include "bonded_pair_particle_arbitrary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "triplet_calc_angular_f.h"


#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()


BondedPairParticleArbitrary::BondedPairParticleArbitrary(Simulation* parent)
  : BondedPairParticleCalc(parent)
{
  init();
}

BondedPairParticleArbitrary::~BondedPairParticleArbitrary()
{
}

void BondedPairParticleArbitrary::init()
{
//   m_properties.setClassName("BondedPairParticleArbitrary");
  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("BondedPairParticleArbitrary");

//   STRINGPC(listName, m_listName, "Identifier of the list of bonded pairs, this Calculator should belong to.");

//   STRINGPC
//       (symbol, m_symbolName,
//        "Name of the symbol for the calculated property.");

  STRINGPC
      (expression, m_expression,
       "The mathematical expression to be computed for the pairFactor."
      );
  
  STRINGPC
      (particleFactor_i, m_1stPExpression,
       "The mathematical expression of the additional particle factor for the first particle."
      );

  STRINGPC
      (particleFactor_j, m_2ndPExpression,
       "The mathematical expression of the additional particle factor for the second particle."
      );

  STRINGPC
      (useOldFor, m_oldSymbols,
       "Here, you can list the used symbols, which should be treated as \"old\", i.e., this calculator will not wait for those symbols to be computed beforehand, but it will take what it finds. Separate the symbols by the \"|\"- (\"pipe\"-) symbol."
      );
    

  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";

// #ifdef _OPENMP
//   m_particleCalculator = true;
// #endif
  m_expression = "idVec(1)";

//   m_listName = "undefined";

  m_oldSymbols = "---";

}

void BondedPairParticleArbitrary::setup()
{
  if(m_expression == "undefined")
    throw gError("BondedPairParticleArbitrary::setup", ": Attribute 'expression' has value \"undefined\"");
       
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  
  m_function.setExpression(m_expression);
  m_function.setColourPair(cp);
  m_1stparticleFactor.setExpression(m_1stPExpression);
  m_1stparticleFactor.setColourPair(cp);
  m_2ndparticleFactor.setExpression(m_2ndPExpression);
  m_2ndparticleFactor.setColourPair(cp);

  BondedPairParticleCalc::setup();
  
}

void BondedPairParticleArbitrary::copyMembersTo(ValCalculator* vc)
{
  BondedPairParticleCalc::copyMembersTo(vc);

   ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  ((BondedPairParticleArbitrary*) vc)->m_function.setExpression(m_expression);
  ((BondedPairParticleArbitrary*) vc)->m_function.setColourPair(cp);
  //       ((BondedPairParticleArbitrary*) vc)->m_function.setReturnType(funcType);
  ((BondedPairParticleArbitrary*) vc)->m_1stparticleFactor.setExpression(m_1stPExpression);
  ((BondedPairParticleArbitrary*) vc)->m_1stparticleFactor.setColourPair(cp);
  ((BondedPairParticleArbitrary*) vc)->m_2ndparticleFactor.setExpression(m_2ndPExpression);
  ((BondedPairParticleArbitrary*) vc)->m_2ndparticleFactor.setColourPair(cp);

}

bool BondedPairParticleArbitrary::findStage()
{
  if(m_stage == -1)
  {
    // no symbols used?
    if(m_function.usedSymbols().empty() && m_1stparticleFactor.usedSymbols().empty() && m_2ndparticleFactor.usedSymbols().empty())
    {
      m_stage = 0;
      MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': no symbols used: stage is now " << m_stage);
      return true;
    }
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;
         
    // loop over the used symbols
    
    typed_value_list_t usedSymbols;
    for(typed_value_list_t::const_iterator s = m_function.usedSymbols().begin(); s != m_function.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    for(typed_value_list_t::const_iterator s = m_1stparticleFactor.usedSymbols().begin(); s != m_1stparticleFactor.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    for(typed_value_list_t::const_iterator s = m_2ndparticleFactor.usedSymbols().begin(); s != m_2ndparticleFactor.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    
    // go through the string of "old" symbols and remove those from usedSymbols 
    if (m_oldSymbols != "---") 
    {
      bool run = true;
      string working = m_oldSymbols;
      while(run)
      {
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
        for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end(); ++s)
        {
          if((*s)->name() == cur)
            toRemove.push_back(*s);
        }
        // remove it
        for(typed_value_list_t::const_iterator s = toRemove.begin(); s != toRemove.end(); ++s)
          usedSymbols.remove(*s);
        if(toRemove.empty())
          throw gError("BondedPairParticleArbitrary::findStage", "Unable to find old symbol '" + cur + "' among used symbols");
      }
    }

    
    for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
    {
      
      string name = (*s)->name();
      MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << ": now checking for used symbol with complete name " << name);
      // unfortunately, we have to remove the brackets for vectors and tensors
      if(name[0] == '{' || name[0] == '[')
      {
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
//       MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << ":short name: " << name);
      
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,
        
//       MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);

      assert(m_phaseUser == 1 || m_phaseUser == 2);

      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
	(
	 M_MANAGER, 
	 vector<ValCalculator*>* vCs;

	 // loop over non-bonded ValCalculators
	 vCs = &(cp->valCalculatorsFlat());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("BondedPairParticleArbitrary::findStage", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	   } // end of loop over non-bonded vCs of current ColourPair

	 // loop over bonded ValCalculators
	 vCs = &(cp->bondedValCalculatorsFlat());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("BondedPairParticleArbitrary::findStage", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	   } // end of loop over bonded vCs of current ColourPair
	 
	 // can we stop FOR_EACH_COLOUR_PAIR now?
	 if(tooEarly) {
	   __cp = __end;
	   // important because there still comes the ++__cp from the loop
	   --__cp;
	 }
	 );
        
          if(!tooEarly)
          {
          // we have to search now in the ParticleCaches
            vector<ParticleCache*>* pCs;
            pCs = &(Particle::s_cached_flat_properties);
            FOR_EACH
	      (vector<ParticleCache*>,
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
			   MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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

          } // end if(!tooEarly) (for search in ParticleCaches)
        
          if(!tooEarly) {
	    // we have to search now in the TripletCalculators
	    vector<TripletCalculator*>* tCs;
	    tCs = M_PHASE->bondedTripletCalculatorsFlat();
	    // MSG_DEBUG("BondedPairParticleArbitrary::findStage", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<TripletCalculator*>, (*tCs),              
	       list<string> symbols = (*__iFE)->mySymbolNames();            
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
		     tooEarly = true;
		     m_stage = -1;
		   }
		   else {
		     if(stage >= m_stage) m_stage = stage+1;
		   }
		 }
	       }
	       // may we abort the loop over the triplet Calculators?
	       if(tooEarly) {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in tripletCalculators)

    } // end loop over used symbols
    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true;
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}


bool BondedPairParticleArbitrary::findStage_0()
{
  if(m_stage == -1)
  {
    // no symbols used?
    if(m_function.usedSymbols().empty() && m_1stparticleFactor.usedSymbols().empty() && m_2ndparticleFactor.usedSymbols().empty())
    {
      m_stage = 0;
      MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': no symbols used: stage is now " << m_stage);
      return true;
    }
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;
         
    // loop over the used symbols
    
    typed_value_list_t usedSymbols;
    for(typed_value_list_t::const_iterator s = m_function.usedSymbols().begin(); s != m_function.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    for(typed_value_list_t::const_iterator s = m_1stparticleFactor.usedSymbols().begin(); s != m_1stparticleFactor.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    for(typed_value_list_t::const_iterator s = m_2ndparticleFactor.usedSymbols().begin(); s != m_2ndparticleFactor.usedSymbols().end() && !tooEarly; ++s)
      usedSymbols.push_back(*s);
    
    // go through the string of "old" symbols and remove those from usedSymbols 
    if (m_oldSymbols != "---") 
    {
      bool run = true;
      string working = m_oldSymbols;
      while(run)
      {
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
        for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end(); ++s)
        {
          if((*s)->name() == cur)
            toRemove.push_back(*s);
        }
        // remove it
        for(typed_value_list_t::const_iterator s = toRemove.begin(); s != toRemove.end(); ++s)
          usedSymbols.remove(*s);
        if(toRemove.empty())
          throw gError("BondedPairParticleArbitrary::findStage_0", "Unable to find old symbol '" + cur + "' among used symbols");
      }
    }

    
    for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
    {
      
      string name = (*s)->name();
      MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << ": now checking for used symbol with complete name " << name);
      // unfortunately, we have to remove the brackets for vectors and tensors
      if(name[0] == '{' || name[0] == '[')
      {
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
//       MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << ":short name: " << name);
      
        // we do if((*vCIt) != this) now, so, next should be unnecessary
/*      if(m_symbolName == name)
      throw gError("BondedPairParticleArbitrary::findStage_0", className() + " reporting: I cannot use my own symbol as argument.");*/
          
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
	 vCs = &(cp->valCalculatorsFlat_0());
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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


	 // loop over bonded ValCalculators
	 vCs = &(cp->bondedValCalculatorsFlat_0());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("BondedPairParticleArbitrary::findStage", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairParticleArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	   } // end of loop over bonded vCs of current ColourPair
	 
	 // can we stop FOR_EACH_COLOUR_PAIR now?
	 if(tooEarly)
	   {
	     __cp = __end;
	     // important because there still comes the ++__cp from the loop
	     --__cp;
	   }
	 );
        
          if(!tooEarly)
          {
          // we have to search now in the ParticleCaches
            vector<ParticleCache*>* pCs;
            pCs = &(Particle::s_cached_flat_properties_0);
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
			   MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
          } // end if(!tooEarly) (for loop over ParticleCaches)

          if(!tooEarly) {
	    // we have to search now in the TripletCalculators
	    vector<TripletCalculator*>* tCs;
	    tCs = M_PHASE->bondedTripletCalculatorsFlat_0();
	    // MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<TripletCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
		     tooEarly = true;
		     m_stage = -1;
		   }
		   else {
		     if(stage >= m_stage) m_stage = stage+1;
		   }
		 }
	       }
	       // may we abort the loop over the triplet Calculators?
	       if(tooEarly) {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in tripletCalculators)

    } // end loop over symbols
    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true;
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("BondedPairParticleArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}
