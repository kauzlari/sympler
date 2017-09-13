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



#include "val_calculator_arbitrary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"


// const SymbolRegister<ValCalculatorRho> val_calc_rho("ValCalculatorRho");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

ValCalculatorArbitrary::ValCalculatorArbitrary(/*Node*/Simulation* parent)
  : NonBondedPairParticleCalculator(parent)
{
  m_datatype = DataFormat::DOUBLE;
/*  m_1stparticleFactor.setReturnType(Variant::SCALAR);
  m_2ndparticleFactor.setReturnType(Variant::SCALAR);*/
  init();
}    

void ValCalculatorArbitrary::init()
{
  m_properties.setClassName("ValCalculatorArbitrary");

  m_properties.setDescription("Computes completely user-defined properties for the particles, which need pair summation. The definition is done with the attribute 'expression'.");
  
  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the calculated property.");

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
    
  INTPC
      (symmetry, m_symmetry, -2,
       "How does the pair expression behave under index interchange? \"1\" means symmetric behaviour. \"-1\" means antisymmetric behaviour.");

  BOOLPC
      (overwrite, m_overwrite,
       "Is this calculator allowed to overwrite already existing symbols " 
           "with name 'symbol' ?");

  DOUBLEPC
      (cutoff, m_cutoff, 0,
       "Specifies the range of the pair summation.");

    
  m_expression = "undefined";
  m_symbolName = "undefined";
  m_oldSymbols = "---";
  m_symmetry = 1;
  m_overwrite = false;
  m_cutoff = 0;

}

void ValCalculatorArbitrary::setup()
{
  MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": START, cutoff = " << m_cutoff);    
    
  if(m_cutoff <= 0)
    throw gError("ValCalculatorArbitrary::setup", className() + ": Attribute 'cutoff' has no value > 0!");
     
  if(m_expression == "undefined")
    throw gError("ValCalculatorArbitrary::setup", className() + ": Attribute 'expression' has value \"undefined\"");
     
  if(m_symmetry != -1 && m_symmetry != 1)
    throw gError("ValCalculatorArbitrary::setup", className() + ": Attribute 'symmetry' may only be \"-1\" (antisymmetry) or \"1\" (symmetry).");
  
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("ValCalculatorArbitrary::setup", className() + ": Attribute 'stage' is " + ObjToString(m_phaseUser) + ", which is none of the allowed values \"0\", \"1\", \"2\".");
   
  if(m_allPairs)
  {
    
    size_t colour;
    
    if(!m_overwrite)
    {
      // is the symbol already existing somewhere?
      for (colour = 0; colour < M_MANAGER->nColours(); ++colour)
      {
        if(Particle::s_tag_format[colour].attrExists(m_symbolName))
          throw gError("ValCalculatorArbitrary::setup", className() + ": Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(colour) + "'. Second definition is not allowed for overwrite = \"no\"");
      }
    }
    colour = 0;

    if(m_overwrite)
    {
      try
      {
        m_slots.first = Particle::s_tag_format[colour].indexOf(m_symbolName, m_datatype);
        m_slots.first = Particle::s_tag_format[colour].offsetByIndex(m_slots.first);
      }
      catch(gError& err)
      {
        throw gError("ValCalculatorArbitrary::setup", className() + ": search for symbol for species '" + M_MANAGER->species(colour) + " failed. The message was " + err.message()); 
      }
    }
    // see CONVENTION5 for rule about persistencies
    else // m_overwrite = false
    {
      m_slots.first = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
    }
    m_slots.second = m_slots.first;
    
    FOR_EACH_COLOUR_PAIR
        (
        M_MANAGER,
      
    // see CONVENTION5 for rule about persistencies
    
    cp->setCutoff(m_cutoff);
    cp->setNeedPairs(true);
    
    if(colour < cp->firstColour() || colour < cp->secondColour())
    {
      ++colour;
      if(cp->firstColour() == colour)
      {
        if(m_overwrite)
        {
          try
          {
            m_slots.first = Particle::s_tag_format[cp->firstColour()].indexOf(m_symbolName, m_datatype);
            m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByIndex(m_slots.first);
          }
          catch(gError& err)
          {
            throw gError("ValCalculatorArbitrary::setup", className() + ": search for symbol for species '" + M_MANAGER->species(cp->firstColour()) + " failed. The message was " + err.message()); 
          }
        }
        else
          // see CONVENTION5 for rule about persistencies
          m_slots.first = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
      }
      if(cp->secondColour() == colour)
      {
        if(m_overwrite)
        {
          try
          {
            m_slots.second = Particle::s_tag_format[cp->secondColour()].indexOf(m_symbolName, m_datatype);
            m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByIndex(m_slots.second);
          }
          catch(gError& err)
          {
            throw gError("ValCalculatorArbitrary::setup", className() + ": search for symbol for species '" + M_MANAGER->species(cp->secondColour()) + " failed. The message was " + err.message()); 
          }
	}
        else
          // see CONVENTION5 for rule about persistencies
          m_slots.second = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persist.second*/false, m_symbolName).offset;
      }
    }
    else
    {
    // here we have to make sure at least that the offsets are correct
      m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByName(m_symbolName);
      m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByName(m_symbolName);
    }
    
    vector<ColourPair*>::iterator cpTester = __cp;
    // is it the last calculator to be created?
    if(++cpTester == __end)
    {
      m_function.setExpression(m_expression);
      m_function.setColourPair(cp);
//       m_function.setReturnType(funcType);
      m_1stparticleFactor.setExpression(m_1stPExpression);
      m_1stparticleFactor.setColourPair(cp);
      m_2ndparticleFactor.setExpression(m_2ndPExpression);
      m_2ndparticleFactor.setColourPair(cp);
/*      m_1stparticleFactor.setReturnType(pFType);
      m_2ndparticleFactor.setReturnType(pFType);*/
      MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": registering me, CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");
      
      if(m_phaseUser == 0)    
        cp->registerCalc_0(this);
      else if(m_phaseUser == 1)    
        cp->registerCalc(this);
      else // so it is 2
      {
        ValCalculator* vc = copyMySelf() /*new ValCalculatorArbitrary(*this)*/;
        assert(((ValCalculatorArbitrary*) vc)->m_symbolName == m_symbolName);
        assert(vc->stage() == m_stage);
        assert(((ValCalculatorArbitrary*) vc)->m_slots.first == m_slots.first);
        assert(((ValCalculatorArbitrary*) vc)->m_slots.second == m_slots.second);
        assert(((ValCalculatorArbitrary*) vc)->m_parent == m_parent);
        assert(((ValCalculatorArbitrary*) vc)->m_overwrite == m_overwrite);
        assert(((ValCalculatorArbitrary*) vc)->m_species.first == m_species.first);
        assert(((ValCalculatorArbitrary*) vc)->m_species.second == m_species.second);
  
        ((ValCalculatorArbitrary*) vc)->m_function.setExpression(m_expression);
        ((ValCalculatorArbitrary*) vc)->m_function.setColourPair(cp);
        ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setExpression(m_1stPExpression);
        ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setColourPair(cp);
        ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setExpression(m_2ndPExpression);
        ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setColourPair(cp);
        assert(((ValCalculatorArbitrary*) vc)->m_function.returnType() == m_function.returnType()); 
      
        MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": registering copy, CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");    
      
        cp->registerCalc(vc);
         
        cp->registerCalc_0(this);
      }
    }
    
    // No? Then make a copy
    else 
    {
      ValCalculator* vc = copyMySelf() /*new ValCalculatorArbitrary(*this)*/;
      assert(((ValCalculatorArbitrary*) vc)->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((ValCalculatorArbitrary*) vc)->m_slots.first == m_slots.first);
      assert(((ValCalculatorArbitrary*) vc)->m_slots.second == m_slots.second);
      assert(((ValCalculatorArbitrary*) vc)->m_parent == m_parent);
      assert(((ValCalculatorArbitrary*) vc)->m_overwrite == m_overwrite);
      assert(((ValCalculatorArbitrary*) vc)->m_species.first == m_species.first);
      assert(((ValCalculatorArbitrary*) vc)->m_species.second == m_species.second);
  
      ((ValCalculatorArbitrary*) vc)->m_function.setExpression(m_expression);
      ((ValCalculatorArbitrary*) vc)->m_function.setColourPair(cp);
//       ((ValCalculatorArbitrary*) vc)->m_function.setReturnType(funcType);
      ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setExpression(m_1stPExpression);
      ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setColourPair(cp);
      ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setExpression(m_2ndPExpression);
      ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setColourPair(cp);
      assert(((ValCalculatorArbitrary*) vc)->m_function.returnType() == m_function.returnType()); 
      
      MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": registering copy, CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");    
      
      if(m_phaseUser == 0)
        cp->registerCalc_0(vc);
      else if(m_phaseUser == 1)
        cp->registerCalc(vc);
      else // so it is 2
      {
        cp->registerCalc(vc);
         
        vc = copyMySelf() /*new ValCalculatorArbitrary(*this)*/;
        assert(((ValCalculatorArbitrary*) vc)->m_symbolName == m_symbolName);
        assert(vc->stage() == m_stage);
        assert(((ValCalculatorArbitrary*) vc)->m_slots.first == m_slots.first);
        assert(((ValCalculatorArbitrary*) vc)->m_slots.second == m_slots.second);
        assert(((ValCalculatorArbitrary*) vc)->m_parent == m_parent);
        assert(((ValCalculatorArbitrary*) vc)->m_overwrite == m_overwrite);
        assert(((ValCalculatorArbitrary*) vc)->m_species.first == m_species.first);
        assert(((ValCalculatorArbitrary*) vc)->m_species.second == m_species.second);
  
        ((ValCalculatorArbitrary*) vc)->m_function.setExpression(m_expression);
        ((ValCalculatorArbitrary*) vc)->m_function.setColourPair(cp);
//       ((ValCalculatorArbitrary*) vc)->m_function.setReturnType(funcType);
        ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setExpression(m_1stPExpression);
        ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setColourPair(cp);
        ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setExpression(m_2ndPExpression);
        ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setColourPair(cp);
        assert(((ValCalculatorArbitrary*) vc)->m_function.returnType() == m_function.returnType()); 
      
        MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": registering copy, CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");    
        
        cp->registerCalc_0(vc);
      }
      
    }
        );
  }
  else /*m_allPairs == false*/
  {
    if(m_species.first == "undefined")
      throw gError("ValCalculatorArbitrary::setup", className() + ": Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled."); 
    if(m_species.second == "undefined")
      throw gError("ValCalculatorArbitrary::setup", className() + ": Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled."); 

    ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);
    
    cp->setCutoff(m_cutoff);
    cp->setNeedPairs(true);

    if(m_overwrite)
    {
      try
      {
        m_slots.first = Particle::s_tag_format[cp->firstColour()].indexOf(m_symbolName, m_datatype);
        m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByIndex(m_slots.first);
      }
      catch(gError& err)
      {
        throw gError("ValCalculatorArbitrary::setup", className() + ": search for symbol for species '" + M_MANAGER->species(cp->firstColour()) + " failed. The message was " + err.message()); 
      }
      try
      {
        m_slots.second = Particle::s_tag_format[cp->secondColour()].indexOf(m_symbolName, m_datatype);
        m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByIndex(m_slots.second);
      }
      catch(gError& err)
      {
        throw gError("ValCalculatorArbitrary::setup", className() + ": search for symbol for species '" + M_MANAGER->species(cp->secondColour()) + " failed. The message was " + err.message()); 
      }      
    }
    else // m_overwrite = false
    {
      if(Particle::s_tag_format[cp->firstColour()].attrExists(m_symbolName)) 
        throw gError("ValCalculatorArbitrary::setup", className() + ": Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(cp->firstColour()) + "'. Second definition is not allowed for overwrite = \"no\"");

      if(Particle::s_tag_format[cp->secondColour()].attrExists(m_symbolName))
        throw gError("ValCalculatorArbitrary::setup", className() + ": Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(cp->secondColour()) + "'. Second definition is not allowed for overwrite = \"no\"");

        
      
      // see CONVENTION5 for rule about persistencies
      m_slots.first = Particle::s_tag_format[cp->firstColour()].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
  
      if(cp->firstColour() != cp->secondColour())
      {
        // see CONVENTION5 for rule about persistencies
        m_slots.second = Particle::s_tag_format[cp->secondColour()].addAttribute(m_symbolName, m_datatype, /*persist.second*/false, m_symbolName).offset;
      }
      else m_slots.second = m_slots.first;
    }
    
//     MSG_DEBUG("ValCalculatorArbitrary::setup", className() + " for " << m_symbolName << ": m_slots = (" << m_slots.first << ", " << m_slots.second << ")");
    
    
    m_function.setExpression(m_expression);
    m_function.setColourPair(cp);
    m_1stparticleFactor.setExpression(m_1stPExpression);
    m_1stparticleFactor.setColourPair(cp);
    m_2ndparticleFactor.setExpression(m_2ndPExpression);
    m_2ndparticleFactor.setColourPair(cp);
    MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": registering, cutoff = " << m_cutoff);
    
    if(m_phaseUser == 0)    
      cp->registerCalc_0(this);
    else if(m_phaseUser == 1)    
      cp->registerCalc(this);
    else // so it is 2
    {
      ValCalculator* vc = copyMySelf() /*new ValCalculatorArbitrary(*this)*/;
      assert(((ValCalculatorArbitrary*) vc)->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((ValCalculatorArbitrary*) vc)->m_slots.first == m_slots.first);
      assert(((ValCalculatorArbitrary*) vc)->m_slots.second == m_slots.second);
      assert(((ValCalculatorArbitrary*) vc)->m_parent == m_parent);
      assert(((ValCalculatorArbitrary*) vc)->m_overwrite == m_overwrite);
      assert(((ValCalculatorArbitrary*) vc)->m_species.first == m_species.first);
      assert(((ValCalculatorArbitrary*) vc)->m_species.second == m_species.second);
  
      ((ValCalculatorArbitrary*) vc)->m_function.setExpression(m_expression);
      ((ValCalculatorArbitrary*) vc)->m_function.setColourPair(cp);
//       ((ValCalculatorArbitrary*) vc)->m_function.setReturnType(funcType);
      ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setExpression(m_1stPExpression);
      ((ValCalculatorArbitrary*) vc)->m_1stparticleFactor.setColourPair(cp);
      ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setExpression(m_2ndPExpression);
      ((ValCalculatorArbitrary*) vc)->m_2ndparticleFactor.setColourPair(cp);
      assert(((ValCalculatorArbitrary*) vc)->m_function.returnType() == m_function.returnType()); 
      
      MSG_DEBUG("ValCalculatorArbitrary::setup", className() + ": registering copy, CP = (" << cp->firstColour() << ", " << cp->secondColour() << ")");    
       
      cp->registerCalc(vc);
      cp->registerCalc_0(this);
    }
  }
}

bool ValCalculatorArbitrary::findStage()
{
  MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << " START: stage = " << m_stage);
  if(m_stage == -1)
  {
    // no symbols used? not overwriting?
    bool usingSymbols = !(m_function.usedSymbols().empty() && m_1stparticleFactor.usedSymbols().empty() && m_2ndparticleFactor.usedSymbols().empty());
    if(!usingSymbols && !m_overwrite)
    {
      m_stage = 0;
      MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': no symbols used: stage is now " << m_stage);
      return true;
    }
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;

    if(usingSymbols) {
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
		throw gError("ValCalculatorArbitrary::findStage", "Unable to find old symbol '" + cur + "' among used symbols");
	    }
	}
      
      for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
	{
	  
	  string name = (*s)->name();
	  MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << ": now checking for used symbol with complete name " << name);
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
	  MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << ":short name: " << name);
	  
	  // we do if((*vCIt) != this) now, so, next should be unnecessary
	  /*      if(m_symbolName == name)
		  throw gError("ValCalculatorArbitrary::findStage", className() + " reporting: I cannot use my own symbol as argument.");*/
	  
	  findStageForSymbolName(name, tooEarly, nothing);
	  
	} // end loop over used symbols
      
    } // end of if(usingSymbols) 

    if(m_overwrite) {
      // if we are overwriting we shouldn't be the first ones to write into this symbol,
      // but should do so at a later stage. So we make an additional check for own symbol(s)
      list<string> symbolNames = mySymbolNames();
      
      for(list<string>::iterator myNamesIt = symbolNames.begin(); myNamesIt != symbolNames.end(); ++myNamesIt)
	
      	findStageForSymbolName(*myNamesIt, tooEarly, nothing);
      
    }

    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true;
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("ValCalculatorArbitrary::findStage", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}


void ValCalculatorArbitrary::findStageForSymbolName(string name, bool& tooEarly, bool& nothing)
{
      // in the following loops we check, whether there exists a
      // ValCalculator for the current symbol. If yes, we check the 
      // stage of it and try to set the stage of this ValCalculator 
      // consistently to it,

      
      MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", className() << " for " << mySymbolName() << ": m_phaseUser = " << m_phaseUser);
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
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", "symbols of VC: ");
		 
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
			     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", "symbols of VC: ");
		 
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
			     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
		     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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


bool ValCalculatorArbitrary::findStage_0()
{
  MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " START: stage = " << m_stage);
  if(m_stage == -1)
  {
    // no symbols used? not overwriting?
    bool usingSymbols = !(m_function.usedSymbols().empty() && m_1stparticleFactor.usedSymbols().empty() && m_2ndparticleFactor.usedSymbols().empty());
    if(!usingSymbols && !m_overwrite)
    {
      m_stage = 0;
      MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': no symbols used: stage is now " << m_stage);
      return true;
    }
  // this is for aborting, when there is a Calculator, which has stage = -1 itself
    bool tooEarly = false;
  // this is for setting the stage to '0' when there is no Calculator at all 
  // (but probably Integrators or s.th. like that)
    bool nothing = true;

    if(usingSymbols) {
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
		throw gError("ValCalculatorArbitrary::findStage_0", "Unable to find old symbol '" + cur + "' among used symbols");
	    }
	}
      
      for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end() && !tooEarly; ++s)
	{
	  
	  string name = (*s)->name();
	  MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << ": now checking for used symbol with complete name " << name);
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
	  MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << ":short name: " << name);
	  
	  // we do if((*vCIt) != this) now, so, next should be unnecessary
	  /*      if(m_symbolName == name)
		  throw gError("ValCalculatorArbitrary::findStage_0", className() + " reporting: I cannot use my own symbol as argument.");*/
          
	  findStageForSymbolName_0(name, tooEarly, nothing);
	  
	} // end loop over symbols
      
    } // end of if(usingSymbols) 
    
    if(m_overwrite) {
      // if we are overwriting we shouldn't be the first ones to write into this symbol,
      // but should do so at a later stage. So we make an additional check for own symbol(s)
      list<string> symbolNames = mySymbolNames();
      
      for(list<string>::iterator myNamesIt = symbolNames.begin(); myNamesIt != symbolNames.end(); ++myNamesIt)
	
      	findStageForSymbolName_0(*myNamesIt, tooEarly, nothing);
      
    }

    if(tooEarly)
      return false;
    if(m_stage == -1)
    {
      if(nothing)
      {
        m_stage = 0;
        MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
        return true;
      } 
      else return false;
    }
    else 
    {
      MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage is now " << m_stage);
      return true;
    }
  } // end if(m_stage == -1)
  else 
  {
    MSG_DEBUG("ValCalculatorArbitrary::findStage_0", className() << " for symbol '"  << m_symbolName << "': stage was already " << m_stage);
    return true;
  }
}


void ValCalculatorArbitrary::findStageForSymbolName_0(string name, bool& tooEarly, bool& nothing)
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
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
		 
		 //             MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", "symbols of VC: ");
		 
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
			     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*vCIt)->className());
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
			   MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
	    // MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<TripletCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
	    // MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", "loop over triplet calculators START");
            FOR_EACH
	      (
	       vector<QuintetCalculator*>, (*tCs),
	       list<string> symbols = (*__iFE)->mySymbolNames();
	       for(list<string>::iterator symIt = symbols.begin(); symIt != symbols.end(); ++symIt) {
		 if(*symIt == name) {
		   nothing = false;
		   int stage = (*__iFE)->stage();
		   if(stage == -1) {
		     MSG_DEBUG("ValCalculatorArbitrary::findStageForSymbolName_0", className() << " for symbol '"  << m_symbolName << "': too early because of " << (*__iFE)->className());
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
