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


#include "symbol.h"
#include "node.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"



#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()

using namespace std;

REGISTER_SMART_ENUM
(SymbolFactory,
 "Expression, computed and saved in a symbol, to be used in expression of other modules.");

Symbol::Symbol(string symbol)
  : Node(NULL), m_stage(-1), m_phase(1), m_symbolName(symbol/*"undefined"*/), m_datatype(DataFormat::DOUBLE), m_overwrite(false)
{
}

Symbol::Symbol(/*Node*/Simulation *parent)
  : Node((Node*) parent), m_stage(-1), m_phase(1)
{
  init();
}
 
/* Methods */

void Symbol::init()
{
  m_properties.setClassName("Symbol");

  INTPC
      (stage, m_phaseUser, -1,
       "When during the timestep should this Calculator be called? Current possibilities are: \n\"0\": At the very beginning of the timestep.\n\"1\": Right before the force calculations.\n\"2\": At both instances.");

//   BOOLPC
//       (overwrite, m_overwrite,
//        "Is this calculator allowed to overwrite already existing symbols " 
//            "with name 'symbol' ?");

  m_symbolName = "undefined";
  m_phaseUser = 1;
  // We set this as the default value for m_overwrite even though only selected 
  // sub-classes provide an XML-attribute for this member
  m_overwrite = false;
}


#ifdef _OPENMP
int Symbol::setNumOfDoubles() {
    return DataFormat::getNumOfDoubles(m_datatype);
}
#endif

void Symbol::checkOverwriteForStageFinding(bool& tooEarly, bool& nothing)
{
    if(m_overwrite) {
      // if we are overwriting we shouldn't be the first ones to write into this symbol,
      // but should do so at a later stage. So we make an additional check for own symbol(s)
      list<string> symbolNames = mySymbolNames();
      
      for(list<string>::iterator myNamesIt = symbolNames.begin(); myNamesIt != symbolNames.end(); ++myNamesIt)
	
      	findStageForSymbolName(*myNamesIt, tooEarly, nothing);
      
    }
}


void Symbol::checkOverwriteForStageFinding_0(bool& tooEarly, bool& nothing)
{
    if(m_overwrite) {
      // if we are overwriting we shouldn't be the first ones to write into this symbol,
      // but should do so at a later stage. So we make an additional check for own symbol(s)
      list<string> symbolNames = mySymbolNames();
      
      for(list<string>::iterator myNamesIt = symbolNames.begin(); myNamesIt != symbolNames.end(); ++myNamesIt)
	
      	findStageForSymbolName_0(*myNamesIt, tooEarly, nothing);
      
    }
}


void Symbol::findStageForSymbolName(string name, bool& tooEarly, bool& nothing)
{
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,

      
      MSG_DEBUG("Symbol::findStageForSymbolName", className() << " for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);
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
		 
		 //             MSG_DEBUG("Symbol::findStageForSymbolName", "symbols of VC: ");
		 
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
			     MSG_DEBUG("Symbol::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		 
		 //             MSG_DEBUG("Symbol::findStageForSymbolName", "symbols of VC: ");
		 
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
			     MSG_DEBUG("Symbol::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	      (
	       vector<ParticleCache*>,
	       (*pCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 {
		   if(*symIt == name)
		     {
		       nothing = false;
		       int stage = (*__iFE)->stage();
		       if(stage == -1) 
			 {
			   MSG_DEBUG("ParticleCacheArbitrary::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
		   // important because there still comes the ++i from the loop
		   --__iFE; 
		 }
	       );
	    
          } // end if(!tooEarly) (for search in ParticleCaches)
        
          if(!tooEarly) {
	    // we have to search now in the TripletCalculators
	    vector<TripletCalculator*>* tCs;
	    tCs = M_PHASE->bondedTripletCalculatorsFlat();
            FOR_EACH
	      (
	       vector<TripletCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("Symbol::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
		 // important because there still comes the ++i from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in tripletCalculators)

          if(!tooEarly) {
	    // we have to search now in the QuintetCalculators
	    vector<QuintetCalculator*>* tCs;
	    tCs = M_PHASE->bondedQuintetCalculatorsFlat();
            FOR_EACH
	      (
	       vector<QuintetCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("Symbol::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
		     tooEarly = true;
		     m_stage = -1;
		   }
		   else {
		     if(stage >= m_stage) m_stage = stage+1;
		   }
		 }
	       }
	       // may we abort the loop over the QuintetCalculators?
	       if(tooEarly) {
		 __iFE = __end;
		 // important because there still comes the ++i from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in QuintetCalculators)
	  
}


void Symbol::findStageForSymbolName_0(string name, bool& tooEarly, bool& nothing)
{
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,
      
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
	(
	 M_MANAGER, 
	 assert(m_phaseUser == 2 || m_phaseUser == 0);
	 // loop over non-bonded ValCalculators
	 vector<ValCalculator*>* vCs;
	 vCs = &(cp->valCalculatorsFlat_0());
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("Symbol::findStageForSymbolName_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("Symbol::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	   } // end loop over non-bonded ValCalculators

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
		 
		 //             MSG_DEBUG("Symbol::findStageForSymbolName_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("Symbol::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	 
	 ); // end loop over ColourPairs
        
          if(!tooEarly)
          {
          // we have to search now in the ParticleCaches
            vector<ParticleCache*>* pCs;
            pCs = &(Particle::s_cached_flat_properties_0);
            FOR_EACH
	      (
	       vector<ParticleCache*>,
	       (*pCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 {
		   if(*symIt == name)
		     {
		       nothing = false;
		       int stage = (*__iFE)->stage();
		       if(stage == -1) 
			 {
			   MSG_DEBUG("Symbol::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
	    // MSG_DEBUG("Symbol::findStageForSymbolName_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<TripletCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("Symbol::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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

          if(!tooEarly) {
	    // we have to search now in the QuintetCalculators
	    vector<QuintetCalculator*>* tCs;
	    tCs = M_PHASE->bondedQuintetCalculatorsFlat_0();
	    // MSG_DEBUG("Symbol::findStageForSymbolName_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<QuintetCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("Symbol::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
		     tooEarly = true;
		     m_stage = -1;
		   }
		   else {
		     if(stage >= m_stage) m_stage = stage+1;
		   }
		 }
	       }
	       // may we abort the loop over the QuintetCalculators?
	       if(tooEarly) {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in QuintetCalculators)

}

