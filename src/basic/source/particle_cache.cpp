/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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



#include "particle_cache.h"

#include "simulation.h"
#include "manager_cell.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


ParticleCache::ParticleCache(/*Node*/Simulation* parent)
  : Symbol(parent)
{
  init();
}

ParticleCache::ParticleCache(size_t colour, size_t offset, string symbolName)
  : Symbol(symbolName)
{
  m_colour = colour;
  m_offset = offset;
}

ParticleCache::~ParticleCache()
{
}

void ParticleCache::init()
{
  m_properties.setClassName("ParticleCache");

  STRINGPC
      (species, m_species,
       "Name for the species of the particles, this Symbol is used for. If set to \"ALL\", the Symbol will be used for all registered species.");
  
  m_species = "undefined";
   
}


void ParticleCache::checkOutputSymbolExistence(size_t colour) {

  if(Particle::s_tag_format[m_colour].attrExists(m_symbolName)) {
    if(m_overwrite) {
      m_offset = Particle::s_tag_format[m_colour].offsetByName(m_symbolName);
      DataFormat::attribute_t tempAttr
	= Particle::s_tag_format[m_colour].attrByName(m_symbolName);
      if(m_datatype != tempAttr.datatype)
	throw gError
	  ("ParticleCache::setup for module " + className(),
	   "Symbol '" + m_symbolName + "' already exists, "
	   "but with different datatype '"
	   + tempAttr.datatypeAsString() + "' to be used due to "
	   "your choice of 'overwrite = \"yes\"', instead of your "
	   "desired datatype '"
	   + DataFormat::attribute_t::datatypeAsString(m_datatype)
	   + "'. Aborting.");
    }
    else
      throw gError
	("ParticleCache::setup for module " + className(), "Symbol '"
	 + m_symbolName + "' was already created by other module "
	 "for species " + M_MANAGER->species(m_colour) + ", and you "
	 "have chosen overwrite = 'no'.");
  } // end of if(Particle::s_tag_format[m_colour].attrExists(..))
  else {
    if(m_overwrite)
      throw gError
	("ParticleCache::setup for module " + className(), "You have "
	 "chosen 'overwrite = \"yes\"', but symbol '" + m_symbolName +
	 "' does not yet exist. Aborting.");
    else
      m_offset = Particle::s_tag_format[m_colour].addAttribute
	(m_symbolName, m_datatype, Symbol::s_persistency, m_symbolName).offset;
  } // end else of if(Particle::s_tag_format[m_colour].attrExists(..))
  
}


void ParticleCache::setup()
{
  Symbol::setup();
  
  if(m_species == "undefined")
    throw gError("ParticleCache::setup for module " + className(), "Attribute 'species' has value \"undefined\"");

  // should we create a Cache for the other colours, too?
  if(m_species == "ALL") {
    // YES! Using m_colour in this loop is correct since we want to
    // create one ParticleCache per colour
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour) {

      checkOutputSymbolExistence(m_colour);

      checkInputSymbolExistences(m_colour);      
      // THE FOLLOWING IS A REMINDER HOW SUCH CHECKS IN CHILD CLASSES COULD LOOK LIKE
      // SINCE STRUCTURE IS ALWAYS THE SAME THE SUBCLASS METHODS COULD CALL BACK A FUNCTION HERE WHICH JUST TAKES THE SPECIFICITIES AS ARGUMENTS (since structure of check is often the same)
      // if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName)) {
      //   if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
      //     throw gError("PCacheIAPWSIF97::setup", "Symbol '" + m_temperatureName + "' already exists as a non-scalar.");
      //   else m_temperatureOffset = Particle::s_tag_format[m_colour].offsetByName(m_temperatureName);
      // } else 
      //     throw gError("PCacheIAPWSIF97::setup", "Symbol '" + m_temperatureName + "' does not exist but required by this module.");

      // if(Particle::s_tag_format[m_colour].attrExists(m_densityName)) {
      //   if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_densityName).datatype)
      //     throw gError("PCacheIAPWSIF97::setup", "Symbol '" + m_densityName + "' already exists as a non-scalar.");
      //   else m_densityOffset = Particle::s_tag_format[m_colour].offsetByName(m_densityName);
      // } else 
      //     throw gError("PCacheIAPWSIF97::setup", "Symbol '" + m_densityName + "' does not exist but required by this module.");


      // is it the last cache to be created?
      if(m_colour == M_MANAGER->nColours()-1)
      {
        if(m_phaseUser == 0)
          Particle::registerCache_0(this);
        else if(m_phaseUser == 1) 
          Particle::registerCache(this);
        else // so it is 2
        {
	  // FIXME: check if we should test all the copies below for
	  // correctness (all members correctly copied). If yes, write
	  // unittest(s)
          ParticleCache* pc = copyMySelf();
          Particle::registerCache(pc);
          Particle::registerCache_0(this);
        }
      }
      // No? Then make a copy
      else 
      {
        ParticleCache* pc = copyMySelf();

        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);          
          pc = copyMySelf();
          Particle::registerCache_0(pc);
        }
      }
    } // end: for(m_colour = 0;...)
  } // end: if(m_species == "ALL")
  else { /*it is a Symbol limited to one colour*/
    
    m_colour = M_MANAGER->getColour(m_species);
    
    checkOutputSymbolExistence(m_colour);
    
    checkInputSymbolExistences(m_colour);
    // SEE m_species == "all" case above for example-impl in subclass
   
    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else
    {
      ParticleCache* pc = copyMySelf();
      Particle::registerCache(pc);
      Particle::registerCache_0(this);
    }     
  } // end: else of if(m_species == "ALL") 
}

  
void ParticleCache::cleanSymbol(string& name) const
{
  if(name[0] == '{' || name[0] == '[') {
    // remove the first bracket
    name.erase(0, 1);
    // remove the last bracket; don't know why, but with these arguments it works
    name.erase(name.size()-1, name.size()-1);
  }
  MSG_DEBUG("ValCalculator::cleanPairSymbol", className() << ": shortened name of symbol: " << name);
}
