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


#include "pca_matrix_inverse.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"


const SymbolRegister<PCaMatrixInverse> pca_matrix_inverse("MatrixInverse");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCaMatrixInverse::PCaMatrixInverse(/*Node*/Simulation* parent)
  : ParticleCache(parent)
{
  m_datatype = DataFormat::TENSOR;
  m_working_mat = gsl_matrix_alloc(SPACE_DIMS, SPACE_DIMS);

  m_inverse.size1 = m_inverse.size2 = m_inverse.tda = SPACE_DIMS;
  m_inverse.owner = 0;
  m_inverse.data = NULL;
  
  m_permutation = gsl_permutation_alloc(3);
    
  init(); 
}

PCaMatrixInverse::~PCaMatrixInverse()
{
  gsl_matrix_free(m_working_mat);
  gsl_permutation_free(m_permutation);
}

void PCaMatrixInverse::init()
{
  m_properties.setClassName("PCaMatrixInverse");
  m_properties.setName("MatrixInverse");

  m_properties.setDescription("Computes the inverse of a given symmetric square matrix.");
      
  STRINGPC
      (symbol, m_symbolName,
       "Symbol name for the inverse.");
  
  STRINGPC
      (tensor, m_tensor_symbol,
       "Name of the tensor for which to compute the inverse.");
  
  m_symbolName = "undefined";
  m_tensor_symbol = "undefined";   
}

void PCaMatrixInverse::setup()
{
  // turns of the default gsl_error handling so that we can handle ourselves
  gsl_set_error_handler_off ();

  if(m_species == "undefined")
    throw gError("PCaMatrixInverse::setup", "Attribute 'species' has value \"undefined\"."); 

  if(m_symbolName == "undefined")
    throw gError("PCaMatrixInverse::setup", "Attribute 'symbol' has value \"undefined\"."); 

  if(m_tensor_symbol == "undefined")
    throw gError("PCaMatrixInverse::setup", "Attribute 'tensor' has value \"undefined\"."); 

  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("PCaMatrixInverse::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  
  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours too?
  if(m_species == "ALL")
  {
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      // is the symbol already existing somewhere?
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
        throw gError("PCaMatrixInverse::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed in MatrixInverse.");
      // the tensor MUST already exist
      if(!Particle::s_tag_format[m_colour].attrExists(m_tensor_symbol))
        throw gError("PCaMatrixInverse::setup", "Symbol " + m_tensor_symbol + " not found for species '" + M_MANAGER->species(m_colour) + "'.");
    }
    m_colour = 0;
    if(Particle::s_tag_format.size() > 1)
    {
      // loop over colours, except the last one for making copies
      for(m_colour = 0; m_colour < Particle::s_tag_format.size()-1; ++m_colour)
      {
        // FIXME: persistency currently always false -> is this OK?
        m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
        m_tensor_offset = Particle::s_tag_format[m_colour].offsetByName(m_tensor_symbol);
                
        // register a copy
        ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
        assert(((PCaMatrixInverse*) pc)->mySymbolName() == m_symbolName);
        assert(((PCaMatrixInverse*) pc)->stage() == m_stage);
        assert(((PCaMatrixInverse*) pc)->m_colour == m_colour);
        assert(((PCaMatrixInverse*) pc)->m_offset == m_offset);
        
        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);
          
          pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
          assert(((PCaMatrixInverse*) pc)->mySymbolName() == m_symbolName);
          assert(((PCaMatrixInverse*) pc)->stage() == m_stage);
          assert(((PCaMatrixInverse*) pc)->m_colour == m_colour);
          assert(((PCaMatrixInverse*) pc)->m_offset == m_offset);
           
          Particle::registerCache_0(pc);
        }
      }
      // now, m_colour = last colour 
    }
  }
  else
  {
    m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
    // are the symbols already existing?
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      throw gError("PCaMatrixInverse::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed in MatrixInverse.");
    // the tensor MUST already exist
    if(!Particle::s_tag_format[m_colour].attrExists(m_tensor_symbol))
      throw gError("PCaMatrixInverse::setup", "Symbol " + m_tensor_symbol + " not found for species '" + M_MANAGER->species(m_colour) + "'.");
  }
  // next lines are done in any case
  // FIXME: persistency currently always false -> is this OK?
  m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
  
//   MSG_DEBUG("PCaMatrixInverse::setup", "m_evecs_offset = " << m_evecs_offset);
  
  m_tensor_offset = Particle::s_tag_format[m_colour].offsetByName(m_tensor_symbol);
  
  if(m_phaseUser == 0) 
    Particle::registerCache_0(this);
  else if(m_phaseUser == 1) 
    Particle::registerCache(this);
  else // so it is 2
  {
    // register a copy
    ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
    assert(((PCaMatrixInverse*) pc)->mySymbolName() == m_symbolName);
    assert(((PCaMatrixInverse*) pc)->stage() == m_stage);
    assert(((PCaMatrixInverse*) pc)->m_colour == m_colour);
    assert(((PCaMatrixInverse*) pc)->m_offset == m_offset);
        
    Particle::registerCache(pc);
    Particle::registerCache_0(this);
    
  }
}


bool PCaMatrixInverse::findStage()
{
  return Symbol::findStageNewPrelim();
}


bool PCaMatrixInverse::findStage_0()
{
  return Symbol::findStageNewPrelim_0();
}




// old version to be deleted if new one works
#if 0

bool PCaMatrixInverse::findStage()
{
//   MSG_DEBUG("PCaMatrixInverse::findStage", "START");
  if(m_stage == -1)
  {
  // in the following loops we check, whether there exists a
  // ValCalculator for the tensor symbol. If yes, we check the 
  // stage of it and try to set the stage of this ParticleCache 
  // consistently to it,
    
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;     
   
  // first, loop over ColourPairs
    FOR_EACH_COLOUR_PAIR
      (
       M_MANAGER, 
       assert(m_phaseUser == 2 || m_phaseUser == 1);
       vector<ValCalculator*>* vCs;
       vCs = &(cp->valCalculatorsFlat());
       // loop over ValCalculators
       for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	   (vCIt != vCs->end() && !tooEarly); ++vCIt)
	 {
	   
	   
	   list<string> symbols = (*vCIt)->mySymbolNames();
	   //           MSG_DEBUG("PCaMatrixInverse::findStage", "symbols of VC: ");
	   
	   /*          for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		       cout << *symIt << endl;*/
	   
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if((*symIt) == m_tensor_symbol)
		 {
		   nothing = false;
		   int stage = (*vCIt)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("PCaMatrixInverse::findStage", " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		  if(*symIt == m_tensor_symbol)
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
        // loop over stages
	vector<ParticleCache*>* pCs;
	/*if(m_phaseUser == 1)*/ pCs = &(Particle::s_cached_flat_properties);
	//           if(m_phaseUser == 0) pCs = &(Particle::s_cached_flat_properties_0);
	FOR_EACH
	  (
	   vector<ParticleCache*>,
	   (*pCs),
	   // at the top we have excluded this ParticleCache to enter this loop
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       
	       if ((*symIt) == m_tensor_symbol)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("PCaMatrixInverse::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
      } // end of if(!tooEarly) (for loop over ParticleCaches)
    
    if(!tooEarly)
      {
        // we have to search now in the TripletCalculators
        // loop over stages
	vector<TripletCalculator*>* tCs;
	tCs = M_PHASE->bondedTripletCalculatorsFlat();
	FOR_EACH
	  (
	   vector<TripletCalculator*>,
	   (*tCs),
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if ((*symIt) == m_tensor_symbol)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("PCaMatrixInverse::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
      } // end of if(!tooEarly) (for loop over TripletCalculators)

    if(!tooEarly)
      {
        // we have to search now in the QuintetCalculators
        // loop over stages
	vector<QuintetCalculator*>* tCs;
	tCs = M_PHASE->bondedQuintetCalculatorsFlat();
	FOR_EACH
	  (
	   vector<QuintetCalculator*>,
	   (*tCs),
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if ((*symIt) == m_tensor_symbol)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("PCaMatrixInverse::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
      } // end of if(!tooEarly) (for loop over QuintetCalculators)

    
    if(tooEarly)
      return false;
    if(m_stage == -1)
      {
	if(nothing)
          {
            m_stage = 0;
            MSG_DEBUG("PCaMatrixInverse::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
            return true; 
          }
	else return false;
      }
    else 
      {
	MSG_DEBUG("PCaMatrixInverse::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
	return true;
      }
  } // end if(m_stage == -1)
  else 
    {
      MSG_DEBUG("PCaMatrixInverse::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
      return true;
    }
}

bool PCaMatrixInverse::findStage_0()
{
//   MSG_DEBUG("PCaMatrixInverse::findStage_0", "START");
  if(m_stage == -1)
  {
  // in the following loops we check, whether there exists a
  // ValCalculator for the tensor symbol. If yes, we check the 
  // stage of it and try to set the stage of this ParticleCache 
  // consistently to it,
    
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;     
   
  // first, loop over ColourPairs
    FOR_EACH_COLOUR_PAIR
        (
	 M_MANAGER, 
	 assert(m_phaseUser == 0 || m_phaseUser == 2);
	 vector<ValCalculator*>* vCs;
	 vCs = &(cp->valCalculatorsFlat_0());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     
	     
	     list<string> symbols = (*vCIt)->mySymbolNames();
	     //           MSG_DEBUG("PCaMatrixInverse::findStage_0", "symbols of VC: ");
	     
	     /*          for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
			 cout << *symIt << endl;*/
	     
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 if((*symIt) == m_tensor_symbol)
		   {
		     nothing = false;
		     int stage = (*vCIt)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		  if(*symIt == m_tensor_symbol)
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
        // loop over stages
          vector<ParticleCache*>* pCs;
          pCs = &(Particle::s_cached_flat_properties_0);
	  FOR_EACH
	    (
	     vector<ParticleCache*>,
	     (*pCs),
	     // at the top we have excluded this ParticleCache to enter this loop
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 
		 if ((*symIt) == m_tensor_symbol)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
        }  // end of if(!tooEarly) (for loop over ParticleCaches)
	
    if(!tooEarly)
      {
        // we have to search now in the TripletCalculators
        // loop over stages
	vector<TripletCalculator*>* tCs;
	tCs = M_PHASE->bondedTripletCalculatorsFlat_0();
	FOR_EACH
	  (
	   vector<TripletCalculator*>,
	   (*tCs),
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if ((*symIt) == m_tensor_symbol)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
      } // end of if(!tooEarly) (for loop over TripletCalculators)

    if(!tooEarly)
      {
        // we have to search now in the QuintetCalculators
        // loop over stages
	vector<QuintetCalculator*>* tCs;
	tCs = M_PHASE->bondedQuintetCalculatorsFlat_0();
	FOR_EACH
	  (
	   vector<QuintetCalculator*>,
	   (*tCs),
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if ((*symIt) == m_tensor_symbol)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
      } // end of if(!tooEarly) (for loop over QuintetCalculators)


        if(tooEarly)
          return false;
        if(m_stage == -1)
        {
          if(nothing)
          {
            m_stage = 0;
            MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
            return true; 
          }
          else return false;
        }
        else 
        {
          MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
          return true;
        }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("PCaMatrixInverse::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}
    
#endif // end of #if 0
