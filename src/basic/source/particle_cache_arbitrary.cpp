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
    // so the attribute shouldn't yet exist
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName)) {
      throw gError("ParticleCacheArbitrary::setupOffset: For module " + className(), "Symbol '" + m_symbolName + "' already existing. Second definition is not allowed for 'overwrite = \"no\"'.");
    }
    else {
    // FIXME: let's try if it works for all cases that persistency = false
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
    }
    MSG_DEBUG("ParticleCacheArbitrary::setupOffset: For module " + className(), "Offset for " << m_symbolName << " = " << m_offset);
  }  
}


void ParticleCacheArbitrary::addMyUsedSymbolsTo(typed_value_list_t& usedSymbols) {

  FunctionParser::addToTypedValueList(m_function.usedSymbols(), usedSymbols);      
}


