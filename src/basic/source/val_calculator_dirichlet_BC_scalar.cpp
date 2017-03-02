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


#include "val_calculator_dirichlet_BC_scalar.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"
#include "particle_cache.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"

// #include <utility>

const SymbolRegister<ValCalculatorDirichletBCScalar> val_calc_dirichlet_BC_scalar("DirichletBCScalar");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorDirichletBCScalar::ValCalculatorDirichletBCScalar(string symbol)
  : ValCalculatorBC(symbol)
{
//   MSG_DEBUG("ValCalculatorDirichletBCScalar::ValCalculatorDirichletBCScalar", "CONSTRUCTOR");
}

ValCalculatorDirichletBCScalar::ValCalculatorDirichletBCScalar(/*Node*/Simulation* parent)
  : ValCalculatorBC(parent)
{
  // next already done in upper class
  //  m_stage = 0;
  init();
}

void ValCalculatorDirichletBCScalar::init()
{
  m_properties.setClassName("DirichletBCScalar");

  m_properties.setDescription("Saves the pair-specific value of an arbitrary scalar of the boundary particle used for applying a Dirichlet boundary condition (BC) in each pair of particles. This calculator uses a linear approximation. The actual value of the Dirichlet boundary condition is assumed to be stored in the symbol of the respective boundary particle given by the attribute 'scalar'.");
  
   STRINGPC
       (scalar, m_scalarName,
        "Name of the scalar the BC should be applied to.");
  
  
//   STRINGPC(wallSpecies, m_wallSpecies, 
//            "Species of the wall particles."
//           );

   m_symbolName = "undefined";
   m_scalarName = "undefined";
   m_wallSpecies = "undefined";
  
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorDirichletBCScalar::setup()
{
  if(m_symbolName == "undefined")
    throw gError("ValCalculatorDirichletBCScalar::setup", "Attribute 'symbol' has value \"undefined\" .");

  if(m_scalarName == "undefined")
    throw gError("ValCalculatorDirichletBCScalar::setup", "Attribute 'scalar' has value \"undefined\" .");

  m_datatype = DataFormat::DOUBLE;
  
  ValCalculatorBC::setup();  

  // now (2016-12-14), ValCalculatorBC::setup() has defined m_species
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  DataFormat::attribute_t firstAttr =
    Particle::s_tag_format[cp->firstColour()].attrByName(m_scalarName);
  
  DataFormat::attribute_t secondAttr =
    Particle::s_tag_format[cp->secondColour()].attrByName(m_scalarName);
  
  if(firstAttr.datatype != DataFormat::DOUBLE)
    throw gError("ValCalculatorDirichletBCScalar::setup", "the symbol " + m_scalarName +
		 " is registerd as a non-scalar for species " +
		 cp->manager()->species(cp->firstColour()));
  
  if(secondAttr.datatype != DataFormat::DOUBLE)
    throw gError("ValCalculatorDirichletBCScalar::setup", "the symbol " + m_scalarName +
		 " is registerd as a non-scalar for species " +
		 cp->manager()->species(cp->secondColour()));
     
   m_scalarOffset.first = 
     Particle::s_tag_format[cp->firstColour()].offsetByName/*indexOf*/(m_scalarName);

   m_scalarOffset.second = 
     Particle::s_tag_format[cp->secondColour()].offsetByName/*indexOf*/(m_scalarName);
  
}

// void /*pair<size_t, size_t>*/ ValCalculatorDirichletBCScalar::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
// {
//   MSG_DEBUG("ValCalculatorDirichletBCScalar::setSlot", "CALLED");
//   m_slot = slot = cp->tagFormat().addAttribute
//       ("ValCalculator_" + myName() + "_" + cp->toString(), DataFormat::POINT, false, m_symbolName).offset;
// }


// FIXME: Generalise and extract the findStage methods. It is a mess that every module must implement them
// by itself by mainly copy&paste. At the very least, when introducing also ValCalculatorDirichletBCVector etc.
// You must out this method into a ValCalculatorDirichletBCArbitrary (see e.g., ValCalculatorArbitrary)
bool ValCalculatorDirichletBCScalar::findStage()
{
  MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage","START: stage = " << m_stage);
  if(m_stage == -1)
  {

  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for throwing an exception if it is still true at the end
    bool nothing = true;
         


    // LOOP OVER USED SYMBOLS NOT NECESSARY BECAUSE THIS ValCalculator JUST USES ONE!
//     for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
//     {
      
//      string name = (*s)->name();
//       MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << ": now checking for used symbol with complete name " << name);

          
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,
        
      MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage"," for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);
      assert(m_phaseUser == 1 || m_phaseUser == 2);
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
	(
	 M_MANAGER,
	 // loop over non-bonded ValCalculators
	 vector<ValCalculator*>* vCs;
	 vCs = &(cp->valCalculatorsFlat());
	 // loop over ValCalculators
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", "symbols of VC: ");
		 
		 //             for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 //               cout << *symIt << endl;
		 
                 
		 for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		   {
		     if(*symIt == m_scalarName)
		       {
			 nothing = false;
			 int stage = (*vCIt)->stage();
			 if(stage == -1) 
			   {
			     MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
			     tooEarly = true;
			     m_stage = -1;
			   }
			 else
			   {
			     if(stage >= m_stage) m_stage = stage+1;
			   }
		       }
		   }
	       } // end of if((*vCIt) != this)       
	   } // end of for(vector<ValCalculator*>::iterator...

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
		 
		 //             MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", "symbols of VC: ");
		 
		 //             for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 //               cout << *symIt << endl;
		 
		 
		 for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		   {
		     if(*symIt == m_scalarName /*name*/)
		       {
			 nothing = false;
			 int stage = (*vCIt)->stage();
			 if(stage == -1) 
			   {
			     MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
      
      if(!tooEarly) {
	// we have to search now in the ParticleCaches
	vector<ParticleCache*>* pCs;
	pCs = &(Particle::s_cached_flat_properties);
	FOR_EACH
	  (
	   vector<ParticleCache*>,
	   (*pCs)/*Particle::s_cached_flat_properties*/,
	   
	   list<string> symbols = (*__iFE)->mySymbolNames();
	   
	   //               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", "symbols of Cache: ");
	   /*              for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
			   cout << *symIt << endl;*/
	   
	   for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	     {
	       if(*symIt == m_scalarName)
		 {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) 
		     {
		       MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
		 if(*symIt == m_scalarName /*name*/) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
		 if(*symIt == m_scalarName /*name*/) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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


      
//     } // end loop over symbols WHICH IS NOT NEEDED SINCE THIS VALCALCULATOR JUST USES ONE!

    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        throw gError("ValCalculatorDirichletBCScalar::findStage"," for symbol '"  + m_symbolName + "': Didn't find Calculators computing my input value '" + m_scalarName + "'. Something seems to be messy. If you cannot find a reason then contact the programmers.");
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}

bool ValCalculatorDirichletBCScalar::findStage_0()
{
  MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0","START: stage = " << m_stage);
  if(m_stage == -1)
  {

  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for throwing an exception if it is still true at the end
    bool nothing = true;
         

    
//     for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
//     {
      
//      string name = (*s)->name();
//       MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", className() << ": now checking for used symbol with complete name " << name);

          
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,
        
      MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0"," for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);
      assert(m_phaseUser == 0 || m_phaseUser == 2);
      // first, loop over ColourPairs
      FOR_EACH_COLOUR_PAIR
	(M_MANAGER, 
	 vector<ValCalculator*>* vCs;
	 // loop over non-bonded ValCalculators
	 vCs = &(cp->valCalculatorsFlat_0());
	 for(vector<ValCalculator*>::iterator vCIt = vCs->begin(); 
	     (vCIt != vCs->end() && !tooEarly); ++vCIt)
	   {
	     // we have to exclude this ValCalculator
	     if((*vCIt) != this)
	       {
		 list<string> symbols = (*vCIt)->mySymbolNames();
		 
		 //             MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "symbols of VC: ");
		 
		 //             for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 //               cout << *symIt << endl;
		 
                 
		 for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		   {
		     if(*symIt == m_scalarName)
		       {
			 nothing = false;
			 int stage = (*vCIt)->stage();
			 if(stage == -1) 
			   {
			     MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStage", "symbols of VC: ");
		 
		 //             for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		 //               cout << *symIt << endl;
		 
		 
		 for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
		   {
		     if(*symIt == m_scalarName)
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
	  if(m_phaseUser == 0) pCs = &(Particle::s_cached_flat_properties_0);
	  //             MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "loop over ParticleCaches START");
	  FOR_EACH
	    (
	     vector<ParticleCache*>,
	     (*pCs)/*Particle::s_cached_flat_properties*/,
	     
	     //               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "loop over ParticleCaches. Now at = " << (*i)->mySymbolName() << " for symbol " << name);
	     
             
	     list<string> symbols = (*__iFE)->mySymbolNames();
	     
	     //               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "symbols of Cache: ");
	     
	     /*              for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
			     cout << *symIt << endl;*/
	     
	     for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
	       {
		 if(*symIt == m_scalarName)
		   {
		     nothing = false;
		     int stage = (*__iFE)->stage();
		     if(stage == -1) 
		       {
			 MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
	} // end if(!tooEarly) (for loop over ParticleCaches)


          if(!tooEarly) {
	    // we have to search now in the TripletCalculators
	    vector<TripletCalculator*>* tCs;
	    tCs = M_PHASE->bondedTripletCalculatorsFlat_0();
	    // MSG_DEBUG("ValCalculatorArbitrary::findStage_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<TripletCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == m_scalarName) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
	    // MSG_DEBUG("ValCalculatorArbitrary::findStage_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<QuintetCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == m_scalarName) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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


      //     } // end loop over symbols, which is not necessary because this ValCalculator depends on only one symbol

    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        throw gError("ValCalculatorDirichletBCScalar::findStage_0"," for symbol '"  + m_symbolName + "': Didn't find Calculators computing my input values '" + m_scalarName + "'. Something seems to be messy. If you cannot find a reason then contact the programmers.");
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}
