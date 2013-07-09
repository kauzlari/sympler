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

// #include <utility>

const SymbolRegister<ValCalculatorDirichletBCScalar> val_calc_dirichlet_BC_scalar("ValCalculatorDirichletBCScalar");

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
  m_properties.setClassName("DirichletBCVels");

  m_properties.setDescription("Saves the pair-specific value of an arbitrary scalar of the boundary particle used for applying a Dirichlet boundary condition (BC) in each pair of particles. This calculator uses a linear approximation.");
  
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

  // now, ValCalculatorBC::setup() has defined m_species
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  size_t freeColour;
  string freeSpecies;
  if(m_wallIsSecond)
    {
      freeColour = cp->firstColour();
      freeSpecies = cp->firstSpecies();
    }
  else
    {
      freeColour = cp->secondColour();
      freeSpecies = cp->secondSpecies();
    }
    
    
  if(Particle::s_tag_format[cp->firstColour()].attrExists(m_symbolName))
    throw gError("ValCalculatorDirichletBCScalar::setup", "Symbol " + m_symbolName + " already exists for species '" + cp->firstSpecies() + "'.");
    
  if(Particle::s_tag_format[cp->secondColour()].attrExists(m_symbolName))
    throw gError("ValCalculatorDirichletBCScalar::setup", "Symbol " + m_symbolName + " already exists for species '" + cp->secondSpecies() + "'.");
    
  if(!Particle::s_tag_format[cp->firstColour()].attrExists(m_scalarName))
    throw gError("ValCalculatorDirichletBCScalar::setup", "Symbol " + m_scalarName + " not found for species '" + cp->firstSpecies() + "'.");
    
  if(!Particle::s_tag_format[cp->secondColour()].attrExists(m_scalarName))
    throw gError("ValCalculatorDirichletBCScalar::setup", "Symbol " + m_scalarName + " not found for species '" + cp->secondSpecies() + "'.");
    
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


// Currently inlined
#if 0

#ifndef _OPENMP
    void ValCalculatorDirichletBCScalar::compute(Pairdist* pD)
#else
    void ValCalculatorDirichletBCScalar::compute(Pairdist* pD, size_t thread_no)
#endif
{
  double innerDist;
  double outerDist;

  computeDists(pD, innerDist, outerDist);

    Particle* p1st = pD->firstPart();
    Particle* p2nd = pD->secondPart();

     if(m_wallIsSecond)
     {

           // if all worked fine, we may now compute the value of the outer particle
           // the check for innerDist != HUGE_VAL is done below
       pD->tag.doubleByOffset(m_slot) =
	 (outerDist+innerDist)*(p2nd->tag.doubleByOffset(m_scalarOffset.second))/innerDist-(outerDist/innerDist)*p1st->tag.doubleByOffset(m_scalarOffset.first);

     }
     else // so !m_wallIsSecond
     {

          // if all worked fine, we may now compute the value of the outer particle
          // the check for innerDist != HUGE_VAL is done below
          // the first term is for Dirichlet != 0. We assume the value was assigned 
          // to the wall particles.
        pD->tag.doubleByOffset(m_slot) =
	  (outerDist+innerDist)*(p1st->tag.doubleByOffset(m_scalarOffset.first))/innerDist-(outerDist/innerDist)*p2nd->tag.doubleByOffset(m_scalarOffset.second);

     } // of else of if(m_wallIsSecond)

    if(innerDist == HUGE_VAL)
      throw gError("ValCalculatorDirichletBCScalar::compute", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + "c1=" + ObjToString(pD->firstPart()->c) + ", c2=" + ObjToString(pD->secondPart()->c) + ", r1=" + ObjToString(pD->firstPart()->r) + ", r2=" + ObjToString(pD->secondPart()->r));
         
}
#endif

bool ValCalculatorDirichletBCScalar::findStage()
{
  MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage","START: stage = " << m_stage);
  if(m_stage == -1)
  {

  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for throwing an exception if it is still true at the end
    bool nothing = true;
         

    
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
        }       
      }


	  CHECK FOR BONDED AND TRIPLET VCs IS MISSING;

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
            /*if(m_phaseUser == 1)*/ pCs = &(Particle::s_cached_flat_properties);
//             if(m_phaseUser == 0) pCs = &(Particle::s_cached_flat_properties_0);
            //             MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", "loop over ParticleCaches START");
            FOR_EACH
                (
                vector<ParticleCache*>,
            (*pCs)/*Particle::s_cached_flat_properties*/,
              
//               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", "loop over ParticleCaches. Now at = " << (*i)->mySymbolName() << " for symbol " << name);

              
            list<string> symbols = (*i)->mySymbolNames();
              
//               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", "symbols of Cache: ");
               
/*              for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
            cout << *symIt << endl;*/
            
            for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
            {
              if(*symIt == m_scalarName)
              {
                nothing = false;
                int stage = (*i)->stage();
                if(stage == -1) 
                {
                  MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*i)->className());
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
              i = __end;
                // important because there still comes the ++i from the loop
              --i; 
            }
                );

          }

//     } // end loop over symbols

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
          (
          M_MANAGER, 
      vector<ValCalculator*>* vCs;
      vCs = &(cp->valCalculatorsFlat_0());
      // loop over ValCalculators
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
      }


	  CHECK FOR BONDED AND TRIPLET VCs IS MISSING;


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
             if(m_phaseUser == 0) pCs = &(Particle::s_cached_flat_properties_0);
            //             MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "loop over ParticleCaches START");
            FOR_EACH
                (
                vector<ParticleCache*>,
            (*pCs)/*Particle::s_cached_flat_properties*/,
              
//               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "loop over ParticleCaches. Now at = " << (*i)->mySymbolName() << " for symbol " << name);

              
            list<string> symbols = (*i)->mySymbolNames();
              
//               MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", "symbols of Cache: ");
               
/*              for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
            cout << *symIt << endl;*/
            
            for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt)
            {
              if(*symIt == m_scalarName)
              {
                nothing = false;
                int stage = (*i)->stage();
                if(stage == -1) 
                {
                  MSG_DEBUG("ValCalculatorDirichletBCScalar::findStage_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*i)->className());
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
              i = __end;
                // important because there still comes the ++i from the loop
              --i; 
            }
                );

          }

//     } // end loop over symbols

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
