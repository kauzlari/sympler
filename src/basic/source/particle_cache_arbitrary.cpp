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


#include "particle_cache_arbitrary.h"
#include "manager_cell.h"
#include "simulation.h"
#include "colour_pair.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"


#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


ParticleCacheArbitrary::ParticleCacheArbitrary(/*Node*/Simulation* parent/*size_t colour*/)
/*: m_colour(colour), m_stage(0)*/
  : ParticleCache(parent)
{
  init();
}

ParticleCacheArbitrary::~ParticleCacheArbitrary()
{
}

void ParticleCacheArbitrary::init()
{
  m_properties.setClassName("ParticleCacheArbitrary");

  STRINGPC
      (expression, m_expression,
      "The mathematical expression to be computed."
      );
  
  STRINGPC
      (symbol, m_symbolName,
       "Name for the symbol.");
  
  BOOLPC
      (overwrite, m_overwrite,
       "Is this calculator allowed to overwrite already existing symbols " 
           "with name 'symbol' ?");

  m_overwrite = false;
  m_expression = "undefined";
  m_symbolName = "undefined";
}

void ParticleCacheArbitrary::setup()
{
//   Symbol::setup();

/*  m_function.setExpression(m_expression);
  m_function.setColour(m_colour);*/
  
  if(m_expression == "undefined")
    throw gError("ParticleCache::setup", className() + " reports: Attribute 'expression' has value \"undefined\""); 
  if(m_species == "undefined")
    throw gError("ParticleCache::setup", className() + " reports: Attribute 'species' has value \"undefined\""); 
  if(m_symbolName == "undefined")
    throw gError("ParticleCache::setup", className() + " reports: Attribute 'symbol' has value \"undefined\"");
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("ParticleCache::setup", className() + " reports: Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  // should we create a Cache for the other colours too?
  if(m_species == "ALL")
  {
    if(Particle::s_tag_format.size() > 1)
    {
      // loop over colours, m_colour = 0 will be done afterwards
      for(m_colour = 1; m_colour < Particle::s_tag_format.size(); ++m_colour)
      {
        // m_colour will be used (but not modified) in the following function
        setupOffset();
        
          
        
          // register a copy
          ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
          // Absolutely important because setExpression adds the function 
          // to the "toBeCompiled"-list
          ((ParticleCacheArbitrary*) pc)->m_function.setExpression(m_expression);
          ((ParticleCacheArbitrary*) pc)->m_function.setColour(m_colour);
          
          assert(((ParticleCacheArbitrary*) pc)->m_function.returnType() == m_function.returnType()); 
          assert(((ParticleCacheArbitrary*) pc)->mySymbolName() == m_symbolName);
          assert(((ParticleCacheArbitrary*) pc)->stage() == m_stage);
          assert(((ParticleCacheArbitrary*) pc)->m_colour == m_colour);
          assert(((ParticleCacheArbitrary*) pc)->m_offset == m_offset);
          if(m_phaseUser == 0)
            Particle::registerCache_0(pc);
          else if(m_phaseUser == 1)
            Particle::registerCache(pc);
          else // so, m_phaseUser = 2
          {
            Particle::registerCache_0(pc);
          // register a second copy
            ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
          // Absolutely important because setExpression adds the function 
          // to the "toBeCompiled"-list
            ((ParticleCacheArbitrary*) pc)->m_function.setExpression(m_expression);
            ((ParticleCacheArbitrary*) pc)->m_function.setColour(m_colour);
          
            assert(((ParticleCacheArbitrary*) pc)->m_function.returnType() == m_function.returnType()); 
            assert(((ParticleCacheArbitrary*) pc)->mySymbolName() == m_symbolName);
            assert(((ParticleCacheArbitrary*) pc)->stage() == m_stage);
            assert(((ParticleCacheArbitrary*) pc)->m_colour == m_colour);
            assert(((ParticleCacheArbitrary*) pc)->m_offset == m_offset);
            Particle::registerCache(pc);
          }
      }
    }
    // for the following registration of "this" ParticleCache
    m_colour = 0;
  }
  else m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
  
  // next lines are done in any case 
  setupOffset();
  
  m_function.setExpression(m_expression);
  m_function.setColour(m_colour);
  
  if(m_phaseUser == 0)
    Particle::registerCache_0(this);
  else if(m_phaseUser == 1)
    Particle::registerCache(this);
  else // so, m_phaseUser = 2
  {
    Particle::registerCache_0(this);
          // register a copy
    ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
          // Absolutely important because setExpression adds the function 
          // to the "toBeCompiled"-list
    ((ParticleCacheArbitrary*) pc)->m_function.setExpression(m_expression);
    ((ParticleCacheArbitrary*) pc)->m_function.setColour(m_colour);
          
    assert(((ParticleCacheArbitrary*) pc)->m_function.returnType() == m_function.returnType()); 
    assert(((ParticleCacheArbitrary*) pc)->mySymbolName() == m_symbolName);
    assert(((ParticleCacheArbitrary*) pc)->stage() == m_stage);
    assert(((ParticleCacheArbitrary*) pc)->m_colour == m_colour);
    assert(((ParticleCacheArbitrary*) pc)->m_offset == m_offset);
    Particle::registerCache(pc);
  }

}

void ParticleCacheArbitrary::setupOffset()
{
  if(m_overwrite)
  { 
    // so the attribute should already exist
    try
    {
      m_offset = Particle::s_tag_format[m_colour].indexOf(m_symbolName, m_datatype);
      m_offset = Particle::s_tag_format[m_colour].offsetByIndex(m_offset);
    }
    catch(gError& err)
    {
      throw gError("ParticleCache::setup", "search for symbol failed. The message was " + err.message()); 
    }
  }
  else
  {
    // so the attribute shoudn't yet exist
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      throw gError("ParticleCache::setup", "Symbol " + m_symbolName + " already existing. Second definition is not allowed for 'overwrite = \"no\"'.");
    else 
      // FIXME: let's try if it works for all cases that persistency = false
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
    MSG_DEBUG("ParticleCacheArbitrary::setupOffset", "Offset for " << m_symbolName << " = " << m_offset);
  } 
}

bool ParticleCacheArbitrary::findStage()
{
  if(m_stage == -1)
  {
    // no symbols used?
    if(m_function.usedSymbols().empty())
    {
      m_stage = 0;
      MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;
         
    // loop over the used symbols
    
    for(typed_value_list_t::const_iterator s = m_function.usedSymbols().begin(); s != m_function.usedSymbols().end() && !tooEarly; ++s)
    {
      
      string name = (*s)->name();
      // unfortunately, we have to remove the brackets for vectors and tensors
      if(name[0] == '{' || name[0] == '[')
      {
        // remove the first bracket
        name.erase(0, 1);
        // remove the last bracket; don't know why, but with these arguments it works
        name.erase(name.size()-1, name.size()-1);
//         remove first and last
      }
      
      // below, we make now if((*i) != this), so, next is hopefully unnecessary
//       if(m_symbolName == name)
//         throw gError("ParticleCacheArbitrary::findStage", className() + " reporting: I cannot use my own symbol as argument.");

      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ParticleCache 
      // consistently to it,
        
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
          (
          M_MANAGER, 

      assert(m_phaseUser == 1 || m_phaseUser == 2);
      vector<ValCalculator*>* vCs;
//           if(m_phaseUser == 1)
      vCs = &(cp->valCalculatorsFlat());
/*          if(m_phaseUser == 0)
      vCs = &(cp->valCalculatorsFlat_0());*/
          // loop over ValCalculators
      for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
          (vCIt != vCs->end() && !tooEarly); ++vCIt)
      {
        list<string> symbols = (*vCIt)->mySymbolNames();

//             MSG_DEBUG("ParticleCacheArbitrary::findStage", "symbols of VC: ");
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
              MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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

	  // loop over bonded ValCalculators
	  vCs = &(cp->bondedValCalculatorsFlat());
	  // loop over ValCalculators
	  for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	      (vCIt != vCs->end() && !tooEarly); ++vCIt)
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
			  MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
			  tooEarly = true;
			  m_stage = -1;
			}
		      else
			{
			  if(stage >= m_stage) m_stage = stage+1;
			}
		    }
		}
	    } // end of loop over non-bonded vCs of current ColourPair


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
	  //             if(m_phaseUser == 1)
	  pCs = &(Particle::s_cached_flat_properties);
	  FOR_EACH
	    (
	     vector<ParticleCache*>,
	     (*pCs),
	     // we have to exclude this ParticleCache, otherwise, infinite loop
	     if((*__iFE) != this)
	       {
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
	       }
// 	     else // (*__iFE) = this
// 	       MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << ": excluding myself in loop over ParticleCaches.");
	     );
	  
	} // end of if(!tooEarly) (for loop over ParticleCaches)
        
      if(!tooEarly)
	{
          // we have to search now in the TripletCalculators
	  vector<TripletCalculator*>* tCs;
	  tCs = M_PHASE->bondedTripletCalculatorsFlat();
	  //             MSG_DEBUG("ParticleCacheArbitrary::findStage", "loop over triplet calculators START");
	  FOR_EACH
	    (
	     vector<TripletCalculator*>,
	     (*tCs),
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
// 	     else // (*__iFE) = this
// 	       MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << ": excluding myself in loop over ParticleCaches.");
	     );
	  
	} // end of if(!tooEarly) (for loop over triplet calculators)

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


    } // end loop over symbols
    if(tooEarly)
    {
      m_stage = -1;
      return false;
    }
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true; 
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("ParticleCacheArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}

bool ParticleCacheArbitrary::findStage_0()
{
  if(m_stage == -1)
  {
    // no symbols used?
    if(m_function.usedSymbols().empty())
    {
      m_stage = 0;
      MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;
         
    // loop over the used symbols
    
    for(typed_value_list_t::const_iterator s = m_function.usedSymbols().begin(); s != m_function.usedSymbols().end() && !tooEarly; ++s)
    {
      
      string name = (*s)->name();
      // unfortunately, we have to remove the brackets for vectors and tensors
      if(name[0] == '{' || name[0] == '[')
      {
        // remove the first bracket
        name.erase(0, 1);
        // remove the last bracket; don't know why, but with these arguments it works
        name.erase(name.size()-1, name.size()-1);
//         remove first and last
      }
      
      // below, we make now if((*i) != this), so, next is hopefully unnecessary
//       if(m_symbolName == name)
//         throw gError("ParticleCacheArbitrary::findStage_0", className() + " reporting: I cannot use my own symbol as argument.");

      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ParticleCache 
      // consistently to it,
        
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
          (
          M_MANAGER, 

      assert(m_phaseUser == 0 || m_phaseUser == 2);
      vector<ValCalculator*>* vCs;
//           if(m_phaseUser == 1)
//       vCs = &(cp->valCalculatorsFlat());
//           if(m_phaseUser == 0)
      vCs = &(cp->valCalculatorsFlat_0());
          // loop over ValCalculators
      for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
          (vCIt != vCs->end() && !tooEarly); ++vCIt)
      {
        list<string> symbols = (*vCIt)->mySymbolNames();

//             MSG_DEBUG("ParticleCacheArbitrary::findStage_0", "symbols of VC: ");               
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
              MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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


	  // loop over bonded ValCalculators
	  vCs = &(cp->bondedValCalculatorsFlat_0());
	  for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	      (vCIt != vCs->end() && !tooEarly); ++vCIt)
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
			  MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
			  tooEarly = true;
			  m_stage = -1;
			}
		      else
			{
			  if(stage >= m_stage) m_stage = stage+1;
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
              (*pCs),
              // we have to exclude this ParticleCache, otherwise, infinite loop
              if((*__iFE) != this)
              {
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
              }
//               else // (*__iFE) = this
//                  MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << ": excluding myself in loop over ParticleCaches.");
                );

          } // end of if(!tooEarly) (for loop over ParticleCaches)

      if(!tooEarly)
	{
          // we have to search now in the TripletCalculators
	  vector<TripletCalculator*>* tCs;
	  tCs = M_PHASE->bondedTripletCalculatorsFlat_0();
	  //             MSG_DEBUG("ParticleCacheArbitrary::findStage_0", "loop over triplet calculators START");
	  FOR_EACH
	    (
	     vector<TripletCalculator*>,
	     (*tCs),
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
		 // important because there still comes the ++i from the loop
		 --__iFE; 
	       }
// 	     else // (*__iFE) = this
// 	       MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << ": excluding myself in loop over ParticleCaches.");
	     );
	  
	} // end of if(!tooEarly) (for loop over triplet calculators)


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
    {
      m_stage = -1;
      return false;
    }
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true; 
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("ParticleCacheArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}
