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


  BOOLPC
      (overwrite, m_overwrite,
       "Is this calculator allowed to overwrite already existing symbols " 
           "with name given in attribute 'symbol' ?");

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
  m_overwrite = false;
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
  m_slot = slot = cp->tagFormat().addAttribute
    // OLD-STYLE. No idea why this was done.
    //    ("BondedPairVector_" + cp->toString(), m_datatype, false, "m_scalar").offset;
    // NEW-STYLE (2014-10-31)
    (m_symbolName, m_datatype, false, m_symbolName).offset;
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

bool BondedPairArbitrary::findStage()
{
  if(m_stage == -1)
  {
    // no symbols used?
    if(m_function.usedSymbols().empty())
    {
      m_stage = 0;
      MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': no symbols used: stage is now " << m_stage);
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
          throw gError("BondedPairArbitrary::findStage", "Unable to find old symbol '" + cur + "' among used symbols");
      }
    }

    
    for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
    {
      
      string name = (*s)->name();
      MSG_DEBUG("BondedPairArbitrary::findStage", className() << ": now checking for used symbol with complete name " << name);
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
//       MSG_DEBUG("BondedPairArbitrary::findStage", className() << ":short name: " << name);
      
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,
        
//       MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);

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
		 
		 //             MSG_DEBUG("BondedPairArbitrary::findStage", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		 
		 //             MSG_DEBUG("BondedPairArbitrary::findStage", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	    // MSG_DEBUG("BondedPairArbitrary::findStage", "loop over triplet calculators START");
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

          if(!tooEarly) {
	    // we have to search now in the QuintetCalculators
	    vector<QuintetCalculator*>* tCs;
	    tCs = M_PHASE->bondedQuintetCalculatorsFlat();
	    // MSG_DEBUG("BondedPairArbitrary::findStage", "loop over quintet calculators START");
            FOR_EACH
	      (
	       vector<QuintetCalculator*>, (*tCs),              
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
	       // may we abort the loop over the QuintetCalculators?
	       if(tooEarly) {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in QuintetCalculators)


    } // end loop over used symbols
    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true;
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("BondedPairArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}


bool BondedPairArbitrary::findStage_0()
{
  if(m_stage == -1)
  {
    // no symbols used?
    if(m_function.usedSymbols().empty())
    {
      m_stage = 0;
      MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': no symbols used: stage is now " << m_stage);
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
          throw gError("BondedPairArbitrary::findStage_0", "Unable to find old symbol '" + cur + "' among used symbols");
      }
    }

    
    for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
    {
      
      string name = (*s)->name();
      MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << ": now checking for used symbol with complete name " << name);
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
//       MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << ":short name: " << name);
      
        // we do if((*vCIt) != this) now, so, next should be unnecessary
/*      if(m_symbolName == name)
      throw gError("BondedPairArbitrary::findStage_0", className() + " reporting: I cannot use my own symbol as argument.");*/
          
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
		 
		 //             MSG_DEBUG("BondedPairArbitrary::findStage_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		 
		 //             MSG_DEBUG("BondedPairArbitrary::findStage_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
	    // MSG_DEBUG("BondedPairArbitrary::findStage_0", "loop over triplet calculators START");
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

          if(!tooEarly) {
	    // we have to search now in the QuintetCalculators
	    vector<QuintetCalculator*>* tCs;
	    tCs = M_PHASE->bondedQuintetCalculatorsFlat_0();
	    // MSG_DEBUG("BondedPairArbitrary::findStage_0", "loop over QuintetCalculators START");
            FOR_EACH
	      (
	       vector<QuintetCalculator*>, (*tCs),
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
	       // may we abort the loop over the QuintetCalculators?
	       if(tooEarly) {
		 __iFE = __end;
		 // important because there still comes the ++__iFE from the loop
		 --__iFE; 
	       }
	       );	    
          } // end if(!tooEarly) (for search in QuintetCalculators)

    } // end loop over symbols
    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true;
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("BondedPairArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}




