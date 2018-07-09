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


#include "colour_pair.h"

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"
// OLD-STYLE: preliminary necessary as long as we don't have the complete hierarchy
// #include "bonded_pair_particle_vector.h"
// NEW_STYLE: 2014-10-31
#include "bonded_pair_arbitrary.h"
#include "bonded_pair_particle_calc.h"


#define M_PHASE  m_manager->phase()
#define M_SIMULATION  ((Simulation*) M_PHASE->parent())
#define M_PAIRCREATOR M_PHASE->pairCreator()

// size_t ColourPair::s_maxStage = 0;


ColourPair::ColourPair(const ColourPair& cp)
// next makes sure that all members of the ColourPair and their members know the
// adress of the COPY, i.e., the NEW ColourPair
  : m_tag_format(cp.m_tag_format), m_freePairs(global::n_threads, PairList(this)),
    m_frozenPairs(global::n_threads, PairList(this))   
                      //  : m_freePairs(cp.m_freePairs, this),
                      // m_frozenPairs(cp.m_frozenPairs, this), 
                      // m_tag_format(cp.m_tag_format)
{
  throw gError("ColourPair", "Do not use the copy-constructor. Contact the programmer.");

}


void ColourPair::registerCalc(pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp)
{
  if(vC->parent())
    throw gError("ColourPair::registerCalc(slots, vC, bool)", "The function shouldn't have been called because vC does belong to a node hierarchy. Contact the programmer.");
  
  string name = string(vC->myName());
  /*size_t*/int stage = vC->stage();
            vector<ValCalculator*>::iterator i = m_valCalculators_flat.begin();
//   vector<ValCalculator*>::iterator i = m_valCalculators[stage].begin();
            bool found = false;
// 		for(size_t i = 0; i < m_valCalculators.size(); ++i)
            while(!found && i != m_valCalculators_flat/*[stage]*/.end())
            {
// 			if(m_valCalculators[i] -> myName() == vC -> myName())
              if((*i)->myName() == name && (*i)->stage() == stage)
              {
                delete vC;
                (*i) -> mySlots(theSlots);
                MSG_DEBUG("ColourPair::registerCalc", 
                          "CP" << toString() << ": using existing ValCalculator '" << name << "', slots=" << theSlots.first << ", " << theSlots.second);
                found = true;
              }
              ++i;
            }
            if(!found)
            {
              setCutoff(vC->cutoff());
              m_valCalculators_flat.push_back(vC);
//     m_valCalculators[stage].push_back(vC);
    
    // letting every valCalculator set the slot by itself makes it possible to return 
    // slots from the Pairdist and the particle_t arrays in the same way
    /*return*/ vC -> setSlots(this, theSlots, oneProp);
               MSG_DEBUG("ColourPair::registerCalc", "CP" << toString() << ": creating new ValCalculator '" << vC -> myName() << "', slots=" << theSlots.first << ", " << theSlots.second << ", requested name was " << name);
            }
}

void ColourPair::registerCalc_0(pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp)
{
  if(vC->parent())
    throw gError("ColourPair::registerCalc_0(slots, vC, bool)", "The function shouldn't have been called because vC does belong to a node hierarchy. Contact the programmer.");
  
  string name = string(vC->myName());
  /*size_t*/int stage = vC->stage();
            vector<ValCalculator*>::iterator i = m_valCalculators_flat_0.begin();
//   vector<ValCalculator*>::iterator i = m_valCalculators[stage].begin();
            bool found = false;
// 		for(size_t i = 0; i < m_valCalculators.size(); ++i)
            while(!found && i != m_valCalculators_flat_0/*[stage]*/.end())
            {
// 			if(m_valCalculators[i] -> myName() == vC -> myName())
              if((*i)->myName() == name && (*i)->stage() == stage)
              {
                delete vC;
                (*i) -> mySlots(theSlots);
                MSG_DEBUG("ColourPair::registerCalc_0", 
                          "CP" << toString() << ": using existing ValCalculator '" << name << "', slots=" << theSlots.first << ", " << theSlots.second);
                found = true;
              }
              ++i;
            }
            if(!found)
            {
              setCutoff(vC->cutoff());
              m_valCalculators_flat_0.push_back(vC);
//     m_valCalculators[stage].push_back(vC);
    
    // letting every valCalculator set the slot by itself makes it possible to return 
    // slots from the Pairdist and the particle_t arrays in the same way
    /*return*/ vC -> setSlots(this, theSlots, oneProp);
               MSG_DEBUG("ColourPair::registerCalc_0", "CP" << toString() << ": creating new ValCalculator '" << vC -> myName() << "', slots=" << theSlots.first << ", " << theSlots.second << ", requested name was " << name);
            }
}


void ColourPair::registerBondedCalc(pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp)
{
  if(vC->parent())
    throw gError("ColourPair::registerBondedCalc(slots, vC, bool)", "The function shouldn't have been called because vC does belong to a node hierarchy. Contact the programmer.");
  
  string name = string(vC->myName());
  /*size_t*/int stage = vC->stage();
            vector<ValCalculator*>::iterator i = m_bondedValCalculators_flat.begin();
            bool found = false;
            while(!found && i != m_bondedValCalculators_flat.end())
            {
              if((*i)->myName() == name && (*i)->stage() == stage)
              {
                delete vC;
                (*i) -> mySlots(theSlots);
                MSG_DEBUG("ColourPair::registerBondedCalc", 
                          "CP" << toString() << ": using existing ValCalculator '" << name << "', slots=" << theSlots.first << ", " << theSlots.second);
                found = true;
              }
              ++i;
            }
            if(!found)
            {
              setCutoff(vC->cutoff());
              m_bondedValCalculators_flat.push_back(vC);
    
    // letting every valCalculator set the slot by itself makes it possible to return 
    // slots from the Pairdist and the particle_t arrays in the same way
    /*return*/ vC -> setSlots(this, theSlots, oneProp);
               MSG_DEBUG("ColourPair::registerBondedCalc", "CP" << toString() << ": creating new ValCalculator '" << vC -> myName() << "', slots=" << theSlots.first << ", " << theSlots.second << ", requested name was " << name);
            }
}

void ColourPair::registerBondedCalc_0(/*size_t icl,*/ pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp)
{
  if(vC->parent())
    throw gError("ColourPair::registerBondedCalc_0(slots, vC, bool)", "The function shouldn't have been called because vC does belong to a node hierarchy. Contact the programmer.");
  
  string name = string(vC->myName());
  /*size_t*/int stage = vC->stage();
  vector<ValCalculator*>::iterator i = m_bondedValCalculators_flat_0.begin();
            bool found = false;
            while(!found && i != m_bondedValCalculators_flat_0.end())
            {
              if((*i)->myName() == name && (*i)->stage() == stage)
              {
                delete vC;
                (*i) -> mySlots(theSlots);
                MSG_DEBUG("ColourPair::registerBondedCalc_0", 
                          "CP" << toString() << ": using existing ValCalculator '" << name << "', slots=" << theSlots.first << ", " << theSlots.second);
                found = true;
              }
              ++i;
            }
            if(!found)
            {
              setCutoff(vC->cutoff());
              m_bondedValCalculators_flat_0.push_back(vC);
    
    // letting every valCalculator set the slot by itself makes it possible to return 
    // slots from the Pairdist and the particle_t arrays in the same way
    /*return*/ vC -> setSlots(this, theSlots, oneProp);
               MSG_DEBUG("ColourPair::registerBondedCalc_0", "CP" << toString() << ": creating new ValCalculator '" << vC -> myName() << "', slots=" << theSlots.first << ", " << theSlots.second << ", requested name was " << name);
            }
}

void ColourPair::registerCalc(size_t& slot, ValCalculator* vC, bool oneProp)
{
  if(vC->parent())
    throw gError("ColourPair::registerCalc(slot, vC, bool)", "The function shouldn't have been called because vC does belong to a node hierarchy. Contact the programmer.");
  
  string name = string(vC->myName());
  /*size_t*/int stage = vC->stage();
            vector<ValCalculator*>::iterator i = m_valCalculators_flat.begin();
//   vector<ValCalculator*>::iterator i = m_valCalculators[stage].begin();
            bool found = false;
// 		for(i = m_valCalculators.begin(); i != m_valCalculators.end(); ++i)
// 		for(size_t i = 0; i < m_valCalculators.size(); ++i)
            while(!found && i != m_valCalculators_flat/*[stage]*/.end())
            {
    // the ValCalculators treated by this function know their stage already
              if((*i)->myName() == name && (*i)->stage() == stage)
              {
/*				MSG_DEBUG("ColourPair::registerCalc", 
                "returning existing slot for " << name);*/
                delete vC;
// 				m_valCalculators[i] -> mySlot(slot);
                (*i) -> mySlot(slot);
                MSG_DEBUG("ColourPair::registerCalc", 
                          "CP" << toString() << ": returning existing ValCalculator '" << name << "', slot=" << slot);
                found = true;
              }
              ++i;
            }
            if(!found)
            {
              setCutoff(vC->cutoff());
              m_valCalculators_flat.push_back(vC);
//     m_valCalculators[stage].push_back(vC);
    
    // letting every valCalculator set the slot by itself makes it possible to return 
    // slots from the Pairdist and the particle_t arrays in the same way
              vC -> setSlot(this, slot, oneProp);
              MSG_DEBUG("ColourPair::registerCalc", "CP" << toString() << ": creating new ValCalculator '" 
                  << vC -> myName() << "', slot=" << slot << ", name was " << name);
// 			MSG_DEBUG("ColourPair::registerCalc", "vC->m_slot" << ((ValCalculatorPair*) vC) -> m_slot);
            }  
}
	
void ColourPair::registerCalc_0(size_t& slot, ValCalculator* vC, bool oneProp)
{
  if(vC->parent())
    throw gError("ColourPair::registerCalc_0(slot, vC, bool)", "The function shouldn't have been called because vC does belong to a node hierarchy. Contact the programmer.");
  
  string name = string(vC->myName());
  /*size_t*/int stage = vC->stage();
            vector<ValCalculator*>::iterator i = m_valCalculators_flat_0.begin();
//   vector<ValCalculator*>::iterator i = m_valCalculators[stage].begin();
            bool found = false;
// 		for(i = m_valCalculators.begin(); i != m_valCalculators.end(); ++i)
// 		for(size_t i = 0; i < m_valCalculators.size(); ++i)
            while(!found && i != m_valCalculators_flat_0/*[stage]*/.end())
            {
    // the ValCalculators treated by this function know their stage already
              if((*i)->myName() == name && (*i)->stage() == stage)
              {
/*				MSG_DEBUG("ColourPair::registerCalc", 
                "returning existing slot for " << name);*/
                delete vC;
// 				m_valCalculators[i] -> mySlot(slot);
                (*i) -> mySlot(slot);
                MSG_DEBUG("ColourPair::registerCalc_0", 
                          "CP" << toString() << ": returning existing ValCalculator '" << name << "', slot=" << slot);
                found = true;
              }
              ++i;
            }
            if(!found)
            {
              setCutoff(vC->cutoff());
              m_valCalculators_flat_0.push_back(vC);
//     m_valCalculators[stage].push_back(vC);
    
    // letting every valCalculator set the slot by itself makes it possible to return 
    // slots from the Pairdist and the particle_t arrays in the same way
              vC -> setSlot(this, slot, oneProp);
              MSG_DEBUG("ColourPair::registerCalc_0", "CP" << toString() << ": creating new ValCalculator '" 
                  << vC -> myName() << "', slot=" << slot << ", name was " << name);
// 			MSG_DEBUG("ColourPair::registerCalc", "vC->m_slot" << ((ValCalculatorPair*) vC) -> m_slot);
            }  
}
	
// this function assumes that the Valcalculator does all the registering work itself
void ColourPair::registerCalc(ValCalculator* vC)
{
  // this should be a safe check whether this function should be used
  if(!vC->parent())
    throw gError("ColourPair::registerCalc(ValCalculator* vC)", "The function shouldn't have been called because vC does not belong to a node hierarchy. Contact the programmer.");
  
//   m_valCalculators[vC->stage()].push_back(vC);
  m_valCalculators_flat.push_back(vC);
  MSG_DEBUG("ColourPair::registerCalc", "m_valCalculators_flat.size now = " << m_valCalculators_flat.size() << ", name[last] = " << m_valCalculators_flat[m_valCalculators_flat.size()-1]->className());  
}

// this function assumes that the Valcalculator does all the registering work itself
void ColourPair::registerCalc_0(ValCalculator* vC)
{
  // this should be a safe check whether this function should be used
  if(!vC->parent())
    throw gError("ColourPair::registerCalc_0(ValCalculator* vC)", "The function shouldn't have been called because vC does not belong to a node hierarchy. Contact the programmer.");
  
//   m_valCalculators[vC->stage()].push_back(vC);
  m_valCalculators_flat_0.push_back(vC);
  MSG_DEBUG("ColourPair::registerCalc_0", "m_valCalculators_flat_0.size now = " << m_valCalculators_flat_0.size() << ", name[last] = " << m_valCalculators_flat_0[m_valCalculators_flat_0.size()-1]->className());  
}

	
// this function assumes that the Valcalculator does all the registering work itself
void ColourPair::registerBondedCalc(ValCalculator* vC)
{
  // this should be a safe check whether this function should be used
  if(!vC->parent())
    throw gError("ColourPair::registerBondedCalc(index, ValCalculator* vC)", "The function shouldn't have been called because vC does not belong to a node hierarchy. Contact the programmer.");
  
  m_bondedValCalculators_flat.push_back(vC);

//   MSG_DEBUG("ColourPair::registerBondedCalc", "m_bondedValCalculators_flat.size now = " << m_bondedValCalculators_flat.size() << ", name[last] = " << m_bondedValCalculators_flat[m_bondedValCalculators_flat.size()-1]->className());  
}

// this function assumes that the Valcalculator does all the registering work itself
void ColourPair::registerBondedCalc_0(ValCalculator* vC)
{
  // this should be a safe check whether this function should be used
  if(!vC->parent())
    throw gError("ColourPair::registerBondedCalc_0(index, ValCalculator* vC)", "The function shouldn't have been called because vC does not belong to a node hierarchy. Contact the programmer.");
  
  m_bondedValCalculators_flat_0.push_back(vC);

//   MSG_DEBUG("ColourPair::registerBondedCalc_0", "m_bondedValCalculators_flat_0.size now = " << m_bondedValCalculators_flat_0.size() << ", name[last] = " << m_bondedValCalculators_flat_0[m_bondedValCalculators_flat_0.size()-1]->className());  
}

const string &ColourPair::firstSpecies() const {
  return m_manager->species(m_1stColour);
}

const string &ColourPair::secondSpecies() const {
  return m_manager->species(m_2ndColour);
}


void ColourPair::setCutoff(double c) {
  m_cutoff = max(m_cutoff, c);
  M_SIMULATION->maxCutoff = max(M_SIMULATION->maxCutoff, c);
  MSG_DEBUG("ColourPair::setCutoff", "CP(" << m_1stColour << m_2ndColour << "):now cp-rc=" << m_cutoff << ", now sim-maxrc=" << M_SIMULATION->maxCutoff);
}

bool ColourPair::findStages()
{    
//   MSG_DEBUG("ColourPair::findStages", toString());
  
  bool finished = true;
  for(vector<ValCalculator*>::iterator i = m_valCalculators_flat.begin(); i != m_valCalculators_flat.end(); ++i)
      // must be this order so that findSTage() is definitely executed
    finished = (*i)->findStage() && finished; 

  for(vector<ValCalculator*>::iterator i = m_bondedValCalculators_flat.begin(); i != m_bondedValCalculators_flat.end(); ++i)
    // must be this order so that findSTage() is definitely executed
    finished = (*i)->findStage() && finished; 

  return finished;
}

bool ColourPair::findStages_0()
{    
//   MSG_DEBUG("ColourPair::findStages_0", toString());
  
  bool finished = true;
  for(vector<ValCalculator*>::iterator i = m_valCalculators_flat_0.begin(); i != m_valCalculators_flat_0.end(); ++i)
      // must be this order so that findSTage() is definitely executed
    finished = (*i)->findStage_0() && finished; 

  for(vector<ValCalculator*>::iterator i = m_bondedValCalculators_flat_0.begin(); i != m_bondedValCalculators_flat_0.end(); ++i)
    // must be this order so that findSTage() is definitely executed
    finished = (*i)->findStage_0() && finished; 

  return finished;
}

void ColourPair::sortStages()
{
  m_maxStage = 0;

  for(vector<ValCalculator*>::iterator i = m_valCalculators_flat.begin(); i != m_valCalculators_flat.end(); ++i)
  {
    int stage = (*i)->stage();
    assert(stage > -1);
    if((size_t) stage > /*s_maxStage*/ m_maxStage) /*s_maxStage*/m_maxStage = stage;

    if(m_valCalculators.size() <= (size_t) stage)
    {
      m_valCalculators.resize(stage+1);
    }
#ifdef _OPENMP
    if((m_valCalculatorParts.size() <= (size_t) stage))
    {
      m_valCalculatorParts.resize(stage+1);
    }
#endif

#ifndef _OPENMP
    assert((size_t)stage < m_valCalculators.size());
    m_valCalculators[stage].push_back(*i);
#else

    if ((*i)->particleCalculator()) {
      assert((size_t)stage < m_valCalculatorParts.size());
      m_valCalculatorParts[stage].push_back((ValCalculatorPart*)(*i));  
    }
    else {
      assert((size_t)stage < m_valCalculators.size());
      m_valCalculators[stage].push_back(*i);    
    }
#endif
  } // end of loop over m_valCalculators_flat

  if(!m_valCalculators.empty())
  {
    for(size_t s = 0; s <= m_maxStage /*s_maxStage*/; ++s)
    {
      for(vector<ValCalculator*>::iterator i = m_valCalculators[s].begin(); i != m_valCalculators[s].end(); ++i)
      {
        MSG_DEBUG("ColourPair::sortStages", "found stage " << s << " for " << (*i)->className());
      }
    }
  }
  else m_valCalculators.resize(1);
#ifdef _OPENMP
  if(!m_valCalculatorParts.empty())
  {
    for(size_t s = 0; s <= m_maxStage /*s_maxStage*/; ++s)
    {
      for(vector<ValCalculator*>::iterator i = m_valCalculatorParts[s].begin(); i != m_valCalculatorParts[s].end(); ++i)
      {
        MSG_DEBUG("ColourPair::sortStages", "found stage " << s << " for " << (*i)->className());
      }
    }
  }
  else m_valCalculatorParts.resize(1);
#endif

  // bonded Calculators

  m_maxBondedStage.resize(m_connectedLists.size());
  m_bondedValCalculators.resize(m_connectedLists.size());

  for(size_t icl = 0; icl < m_maxBondedStage.size(); ++icl) 
    m_maxBondedStage[icl] = 0;

  for(vector<ValCalculator*>::iterator i = m_bondedValCalculators_flat.begin(); i != m_bondedValCalculators_flat.end(); ++i)
    {
      // NEW STYLE: 2014-10-31. Note that in the "ValCalculatorPair"-case the correct name might change if you start implementing non-arbitrary Calculators one day. Then you probably need a parent "...Calc"-class as in the "BondedPairParticleCalc"-case.
      string listName;
      if((*i) -> className() == "ValCalculatorPart") {
	listName = ((BondedPairParticleCalc*) (*i))->listName();
      }
      else if((*i) -> className() == "ValCalculatorPair") {
	listName = ((BondedPairArbitrary*) (*i))->listName();
      }
      else throw gError("ColourPair::sortStages", "Can not find stage for bonded calculator '" + (*i) -> name() + "' because don't know how to handle its polymorphic class name '" + (*i) -> className() + "'. Looks like a bug for a programmer.");
      // OLD preliminary cast as long as we don't have the complete hierarchy replace by NEW STYLE above
//       string listName = ((BondedPairParticleVector*) (*i))->listName();

      size_t listIndex = connectedListIndex(listName);

      int stage = (*i)->stage();
      assert(stage > -1);

      if((size_t) stage > m_maxBondedStage[listIndex]) m_maxBondedStage[listIndex] = stage;
      
      if(m_bondedValCalculators[listIndex].size() <= (size_t) stage)
	{
	  m_bondedValCalculators[listIndex].resize(stage+1);
	}
#if 0 // this is when we will have a separation of Calculators
#ifdef _OPENMP
      if((m_valCalculatorParts.size() <= (size_t) stage))
	{
	  m_valCalculatorParts.resize(stage+1);
	}
#endif
#endif
      
#ifndef _OPENMP
      assert((size_t)stage < m_bondedValCalculators[listIndex].size());
      m_bondedValCalculators[listIndex][stage].push_back(*i);
#else
      // The commented out code will be needed when we introduce a separation of Calculators
      // 	if ((*i)->particleCalculator()) {
      // 	  assert((size_t)stage < m_valCalculatorParts.size());
// 	  m_valCalculatorParts[stage].push_back((ValCalculatorPart*)(*i));  
// 	}
// 	else {
      assert((size_t)stage < m_bondedValCalculators[listIndex].size());
      m_bondedValCalculators[listIndex][stage].push_back(*i);    
      // 	}
#endif
    } // end of loop over m_bondedValCalculators_flat

  for(size_t listIndex = 0; listIndex < m_bondedValCalculators.size(); ++listIndex) {
    if(!m_bondedValCalculators[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<ValCalculator*>::iterator i = m_bondedValCalculators[listIndex][s].begin(); i != m_bondedValCalculators[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("ColourPair::sortStages", "found stage " << s << " for " << (*i)->className() );
	      }
	  }
      }
    else m_bondedValCalculators[listIndex].resize(1);
#if 0 // next will be needed when separating the bonded Calculators
#ifdef _OPENMP
    if(!m_bondedalCalculatorParts[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<ValCalculator*>::iterator i = m_bondedValCalculatorParts[listIndex][s].begin(); i != m_bondedValCalculatorParts[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("ColourPair::sortStages", "found stage " << s << " for " << (*i)->className() );
	      }
	  }
      }
    else m_bondedValCalculatorParts[listIndex].resize(1);
#endif
#endif
  }
  MSG_DEBUG("ColourPair::sortStages", "(" << firstColour() << "," << secondColour() << "): highest stage needed =  " << m_maxStage);
  
}


void ColourPair::sortStages_0()
{
  m_maxStage_0 = 0;

  for(vector<ValCalculator*>::iterator i = m_valCalculators_flat_0.begin(); i != m_valCalculators_flat_0.end(); ++i)
  {

    int stage = (*i)->stage();
    assert(stage > -1);
    if((size_t) stage > /*s_maxStage*/ m_maxStage_0) /*s_maxStage*/m_maxStage_0 = stage;
    
    if(m_valCalculators_0.size() <= (size_t) stage)
    {
      m_valCalculators_0.resize(stage+1);
    }
#ifdef _OPENMP
    if((m_valCalculatorParts_0.size() <= (size_t) stage) && ((ValCalculatorPart*)(*i)))
    {
      m_valCalculatorParts_0.resize(stage+1);
    }
#endif
      
#ifndef _OPENMP
    assert((size_t)stage < m_valCalculators_0.size());
    m_valCalculators_0[stage].push_back(*i);
#else
    if ((ValCalculatorPart*)(*i)) {
      assert((size_t)stage < m_valCalculatorParts_0.size());
      m_valCalculatorParts_0[stage].push_back((ValCalculatorPart*)(*i));      
    }
    else {
      assert((size_t)stage < m_valCalculators_0.size());
      m_valCalculators_0[stage].push_back(*i);      
    }
#endif
  }

  if(!m_valCalculators_0.empty())
  {
    for(size_t s = 0; s <= m_maxStage_0 /*s_maxStage*/; ++s)
    {
      for(vector<ValCalculator*>::iterator i = m_valCalculators_0[s].begin(); i != m_valCalculators_0[s].end(); ++i)
      {
        MSG_DEBUG("ColourPair::sortStages_0", "found stage " << s << " for " << (*i)->className());
      }
    }

  }
  else m_valCalculators_0.resize(1);
#ifdef _OPENMP
  if(!m_valCalculatorParts_0.empty())
  {
    for(size_t s = 0; s <= m_maxStage_0 /*s_maxStage*/; ++s)
    {
      for(vector<ValCalculator*>::iterator i = m_valCalculatorParts_0[s].begin(); i != m_valCalculatorParts_0[s].end(); ++i)
      {
        MSG_DEBUG("ColourPair::sortStages_0", "found stage " << s << " for " << (*i)->className());
      }
    }

  }
  else m_valCalculatorParts_0.resize(1);
#endif

  // bonded Calculators
  m_maxBondedStage_0.resize(m_connectedLists.size());
  m_bondedValCalculators_0.resize(m_connectedLists.size());
  for(size_t icl = 0; icl < m_maxBondedStage_0.size(); ++icl) 
    m_maxBondedStage_0[icl] = 0;

  for(vector<ValCalculator*>::iterator i = m_bondedValCalculators_flat_0.begin(); i != m_bondedValCalculators_flat_0.end(); ++i)
    {
      // NEW STYLE: 2014-10-31. Note that in the "BondedPairArbitrary"-case the correct name might change if you start implementing non-arbitrary Calculators one day. Then you probably need a parent "...Calc"-class as in the "BondedPairParticleCalc"-case.
      string listName;
      if((*i) -> name() == "BondedPairParticleCalc") {
	listName = ((BondedPairParticleCalc*) (*i))->listName();
      }
      else if((*i) -> name() == "BondedPairArbitrary") {
	listName = ((BondedPairArbitrary*) (*i))->listName();
      }
      else throw gError("ColourPair::sortStages", "Can not find stage for bonded calculator '" + (*i) -> className() + "' because don't know how to handle its polymorphic class name '" + (*i) -> name() + "'. Looks like a bug for a programmer.");
      // OLD preliminary cast as long as we don't have the complete hierarchy replace by NEW STYLE above
//       string listName = ((BondedPairParticleVector*) (*i))->listName();

      size_t listIndex = connectedListIndex(listName);

      int stage = (*i)->stage();
      assert(stage > -1);

      if((size_t) stage > m_maxBondedStage_0[listIndex]) m_maxBondedStage_0[listIndex] = stage;
      
      if(m_bondedValCalculators_0[listIndex].size() <= (size_t) stage)
	{
	  m_bondedValCalculators_0[listIndex].resize(stage+1);
	}
#if 0 // this is when we will have a separation of Calculators
#ifdef _OPENMP
      if((m_valCalculatorParts_0.size() <= (size_t) stage))
	{
	  m_valCalculatorParts_0.resize(stage+1);
	}
#endif
#endif
      
#ifndef _OPENMP
      assert((size_t)stage < m_bondedValCalculators_0[listIndex].size());
      m_bondedValCalculators_0[listIndex][stage].push_back(*i);
#else
      // The commented out code will be needed when we introduce a separation of Calculators
      // 	if ((*i)->particleCalculator()) {
      // 	  assert((size_t)stage < m_valCalculatorParts_0.size());
      // 	  m_valCalculatorParts_0[stage].push_back((ValCalculatorPart*)(*i));  
      // 	}
      // 	else {
      assert((size_t)stage < m_bondedValCalculators_0[listIndex].size());
      m_bondedValCalculators_0[listIndex][stage].push_back(*i);    
      // 	}
#endif
    } // end of loop over m_bondedValCalculators_flat

  for(size_t listIndex = 0; listIndex < m_bondedValCalculators_0.size(); ++listIndex) {
    if(!m_bondedValCalculators_0[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage_0[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<ValCalculator*>::iterator i = m_bondedValCalculators_0[listIndex][s].begin(); i != m_bondedValCalculators_0[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("ColourPair::sortStages_0", "found stage " << s << " for " << (*i)->className());
	      }
	  }
      }
    else m_bondedValCalculators_0[listIndex].resize(1);
  
#if 0 // next will be needed when separating the bonded Calculators
#ifdef _OPENMP
    if(!m_bondedalCalculatorParts_0[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage_0[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<ValCalculator*>::iterator i = m_bondedValCalculatorParts_0[listIndex][s].begin(); i != m_bondedValCalculatorParts_0[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("ColourPair::sortStages_0", "found stage " << s << " for " << (*i)->className());
	      }
	  }
      }
    else m_bondedValCalculatorParts_0[listIndex].resize(1);
#endif
#endif
  }

  MSG_DEBUG("ColourPair::sortStages_0", "(" << firstColour() << "," << secondColour() << "): highest stage needed =  " << m_maxStage_0);
  
}


void ColourPair::addPairToConnectionAndWrite(Particle *first_p, Particle *second_p, size_t listIndex) {

  if(m_1stColour != m_2ndColour)
    throw gError("ColourPair::addPairToConnectionAndWrite", "Calling this function is not allowed for ColourPairs with two different species like this one because writing into a new file currently just works for two identical species. Species1 = '" + firstSpecies() + "', Species2 = '" + secondSpecies() + "'.");

  point_t boxSize = M_PHASE->boundary()->boundingBox().size();
//   MSG_DEBUG("ColourPair::addPairToConnection", "using box size " << boxSize);
  double size;
  dist_t d;
  d.abs_square = 0;
  /* Write the connection into a file */
  ofstream file_stream;
  string file_name;
  
  file_name = M_SIMULATION->name();
  file_name += "_connector.con";
  
  file_stream.open(file_name.c_str(), ios_base::app);
  file_stream.seekp(0, std::ios_base::end);

  // first connection in file? Then write the line "!!!" first
  if(!Phase::connectionDeclarationFinished) {
    Phase::connectionDeclarationFinished = true;
    file_stream << "!!!" << endl;
  }

  file_stream << "pair "<< connectedListName(listIndex) << " " << first_p -> mySlot << " ";
  if (first_p->isFrozen)
    file_stream << "frozen" << " ";
  else
    file_stream << "free" << " ";

//  EXTEND: add previous list-length if from other colour

  file_stream << second_p -> mySlot << " ";
  if (second_p->isFrozen)
    file_stream << "frozen" << " "<< endl;
  else
    file_stream << "free" << " "<< endl;
  file_stream.close();
  
  /* Take care: The order defines the direction d is pointing to. */
  for (int _i = 0; _i < SPACE_DIMS; _i++) {
    d.cartesian[_i] = first_p -> r[_i] - second_p -> r[_i];
    //MSG_DEBUG("ColourPair::addPairToConnectionAndWrite", "cartesian before" << d.cartesian[_i]);
    // periodic BCs
    size = boxSize[_i];
    if(d.cartesian[_i] > 0.5*size) d.cartesian[_i] -= size; 
    if(d.cartesian[_i] < -0.5*size) d.cartesian[_i] += size; 
    d.abs_square += d.cartesian[_i]*d.cartesian[_i];
  }
  d.abs = sqrt(d.abs_square);
  (m_connectedLists[listIndex])->newPair().set(d, first_p, second_p, true, true);
  //MSG_DEBUG("ColourPair::addPairToConnectionAndWrite", "cartesian after" << d.cartesian);
  //MSG_DEBUG("ColourPair::addPairToConnectionAndWrite", "abs after: d = " << d.abs);
}


void ColourPair::addPairToConnection(Particle *arg_first_p, Particle *arg_second_p, size_t listIndex) {

  Particle* first_p;
  Particle* second_p;

  // The order of the particles must correspond to the order of the colours in this ColourPair
  if(arg_first_p->c == m_1stColour && arg_second_p->c == m_2ndColour) {
    first_p = arg_first_p;
    second_p = arg_second_p;
  }
  else if(arg_first_p->c == m_2ndColour && arg_second_p->c == m_1stColour) {
    first_p = arg_second_p;
    second_p = arg_first_p;
  }
  point_t boxSize = M_PHASE->boundary()->boundingBox().size();
//   MSG_DEBUG("ColourPair::addPairToConnection", "using box size " << boxSize);
  double size;
  dist_t d;
  d.abs_square = 0;
  
  /* Take care: The order defines the direction d is pointing to. */
  for (int _i = 0; _i < SPACE_DIMS; _i++) {
    d.cartesian[_i] = first_p -> r[_i] - second_p -> r[_i];
    //MSG_DEBUG("ColourPair::addPairToConnection", "cartesian before" << d.cartesian[_i]);
    // periodic BCs
    size = boxSize[_i];
    if(d.cartesian[_i] > 0.5*size) d.cartesian[_i] -= size; 
    if(d.cartesian[_i] < -0.5*size) d.cartesian[_i] += size; 
    d.abs_square += d.cartesian[_i]*d.cartesian[_i];
  }
  d.abs = sqrt(d.abs_square);
  (m_connectedLists[listIndex])->newPair().set(d, first_p, second_p, true, true);
  //MSG_DEBUG("ColourPair::addPairToConnection", "cartesian after" << d.cartesian);
  //MSG_DEBUG("ColourPair::addPairToConnection", "abs after: d = " << d.abs);
}


void ColourPair::updateConnectedDistances() {
  point_t boxSize = M_PHASE->boundary()->boundingBox().size();
  //MSG_DEBUG("ColourPair::updateConnectedDistances", "using box size " << boxSize);
  double size;
  
  for (vector<PairList*>::iterator cLit = m_connectedLists.begin(); cLit != m_connectedLists.end(); ++cLit) {
    
    for (Pairdist *i = (*cLit)->first(); i != NULL; i = i->next) {
      i->calculateCartDistance();
      // periodic BCs
      point_t& cartesian = i->m_distance.cartesian;
      //MSG_DEBUG("ColourPair::updateConnectedDistances", "cartesian before: " << i->m_distance.cartesian);
      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
	size = boxSize[dir];
	// Concerning periodic BCs, currently (2009/04/13) this is the only place where, in principle, we have to ask ourselves whether to take > or >= (< or <=). But also currently (2009/04/13) this question is irrelevant since any direction must have at least 2 cells.
	if(cartesian[dir] > 0.5*size) cartesian[dir] -= size; 
	if(cartesian[dir] < -0.5*size) cartesian[dir] += size; 
      }
      // absolute values
      i->m_distance.calcAbs();
      //MSG_DEBUG("ColourPair::updateConnectedDistances", "cartesian after: " << i->m_distance.cartesian);
    }
  }
}

PairList* ColourPair::connectedList(string name)
{
  return m_connectedLists[connectedListIndex(name)];
}

/*!
 * Return the index of the connected list from \a m_connectedLists corresponding to the identifier \a name.
 */
size_t ColourPair::connectedListIndex(string name)
{
  int i = -1;

//   MSG_DEBUG("ColourPair::connectedListIndex", "called for name " << name);
// MSG_DEBUG("ColourPair::connectedListIndex", "name list size = " << m_connectedListNames.size());
// MSG_DEBUG("ColourPair::connectedListIndex", "list size = " << m_connectedLists.size());

  for (size_t ii = 0; ii < m_connectedListNames.size(); ++ii) {
    if(m_connectedListNames[ii] == name) {
      i = ii;
      break;
    }
  }
  if (i == -1) 
    throw gError("ColourPair::connectedListIndex"+FILE_INFO, "No list of connections with name " + name + " found for ColourPair " + toString() + "!");
  else return size_t(i);
}
  
size_t ColourPair::createConnectedListIndexAndWrite(string name)
{
//   MSG_DEBUG("ColourPair::createConnectedListIndexAndWrite", "name: " << name);

  if(!Phase::connectionDeclarationFinished) {
  // search first
  size_t ii;
  for (ii = 0; ii < m_connectedListNames.size(); ++ii) {
    if(m_connectedListNames[ii] == name) {
      return ii;
    }
  }
  // create new
  m_connectedListNames.push_back(name);
  m_connectedLists.push_back(new PairList(this));

  // write in file
  // FIXME: one of at least four places where we create the same file_name from scratch: redundant and error-prone! Replace by member m_connections_file?
  string file_name = M_SIMULATION->name();
  file_name += "_connector.con";
  ofstream file_stream;
  file_stream.open(file_name.c_str(), ios_base::app);
  file_stream.seekp(0, std::ios_base::end);
  file_stream << "pair " << name << " " << firstSpecies() << " " << secondSpecies() << " " << endl;

  return ii;
  }
  else
    throw gError("ColourPair::createConnectedListIndexAndWrite"+FILE_INFO, "Attempt to create a new connected pair list eventhough it seems too late for that. Contact programmer.");
}
  
size_t ColourPair::createConnectedListIndex(string name)
{
//   MSG_DEBUG("ColourPair::createConnectedListIndex", "name: " << name);

  if(!Phase::connectionDeclarationFinished) {
  // search first
  size_t ii;
  for (ii = 0; ii < m_connectedListNames.size(); ++ii) {
    if(m_connectedListNames[ii] == name) {
      return ii;
    }
  }
  // create new
  m_connectedListNames.push_back(name);
  m_connectedLists.push_back(new PairList(this));

  return ii;
  }
  else
    throw gError("ColourPAir::createConnectedListIndex"+FILE_INFO, "Attempt to create a new connected pair list eventhough it seems too late for that. Contact programmer.");
}
  
int ColourPair::searchConnectedListIndex(string name)
{
  for (size_t ii = 0; ii < m_connectedListNames.size(); ++ii) {
    if(m_connectedListNames[ii] == name) return ii;
  }
  return -1;
}
  
  /*!
   * Return the identifier of the connected list from \a m_connectedLists corresponding to the index \a name.
   */
string ColourPair::connectedListName(size_t listIndex)
{
  if(listIndex < m_connectedListNames.size())
    return m_connectedListNames[listIndex];
  else
  throw gError("ColourPair::connectedListName"+FILE_INFO, "requested index " + ObjToString(listIndex)  + " out of bounds.");

}
