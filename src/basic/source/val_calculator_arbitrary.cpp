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



#include "val_calculator_arbitrary.h"
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
  return Symbol::findStageNewPrelim();
}


bool ValCalculatorArbitrary::findStage_0()
{
  return Symbol::findStageNewPrelim_0();
}


void ValCalculatorArbitrary::addMyUsedSymbolsTo(typed_value_list_t& usedSymbols) {

  FunctionParser::addToTypedValueList(m_function.usedSymbols(), usedSymbols);
  FunctionParser::addToTypedValueList(m_1stparticleFactor.usedSymbols(), usedSymbols);
  FunctionParser::addToTypedValueList(/*const typed_value_list_t&*/ m_2ndparticleFactor.usedSymbols(), /*typed_value_list_t&*/ usedSymbols);
      
}
