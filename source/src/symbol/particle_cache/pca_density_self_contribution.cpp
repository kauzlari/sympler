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



#include "pca_density_self_contribution.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<ParticleCacheDensitySelfContribution> pca_density_self("ParticleCacheDensitySelfContribution");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

ParticleCacheDensitySelfContribution::ParticleCacheDensitySelfContribution
  (size_t colour, size_t offset, WeightingFunction *wf, string symbolName)
  : ParticleCache(colour, offset, symbolName), /*m_colour(colour), m_offset(offset), */m_wf(wf)/*, m_symbolName(symbolName)*/
{
  // the next three belong to ParticleCache
/*  m_colour = colour;
  m_offset = offset;
  m_symbolName = symbolName;*/
  
  m_stage = 0;
//   m_datatype = DataFormat::DOUBLE;
}

ParticleCacheDensitySelfContribution::ParticleCacheDensitySelfContribution
    (/*Node*/Simulation* parent)
  : ParticleCache(parent)
{
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  init();
}


ParticleCacheDensitySelfContribution::~ParticleCacheDensitySelfContribution()
{
}

void ParticleCacheDensitySelfContribution::init()
{
  m_properties.setClassName("ParticleCacheDensitySelfContribution");
  m_properties.setName("DensitySC");

  m_properties.setDescription("Contribution of the particle itself to its local density.");
      
  STRINGPC
      (symbol, m_symbolName,
       "Name for the local density.");
  
  STRINGPC
      (weightingFunction, m_weighting_function,
       "Defines the weighting function to be used for the computation of the local density.");
     
  m_weighting_function = "default";
  

}

void ParticleCacheDensitySelfContribution::setup()
{
  if(m_species == "undefined")
    throw gError("ParticleCacheDensitySelfContribution::setup", "Attribute 'species' has value \"undefined\""); 
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("ParticleCacheDensitySelfContribution::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  
  m_wf = M_SIMULATION->findWeightingFunction(m_weighting_function);
  
  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours, too?
  if(m_species == "ALL")
  {
    // This calculator does not care if the symbol is already existing
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      // at least, if the attribute already exists, we preserve the persistency and can hope that another module has set it correctly
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype)
          throw gError("ParticleCacheDensitySelfContribution::setup", "Symbol " + m_symbolName + " already exists as a non-scalar.");
        else m_offset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_symbolName);
      }
      else
      // if it does not yet exist we set the persistency to false
        m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;
    

      // is it the last cache to be created?
      if(m_colour == M_MANAGER->nColours()-1)
      {
        if(m_phaseUser == 0)
          Particle::registerCache_0(this);
        else if(m_phaseUser == 1)
          Particle::registerCache(this);
        else // so it is 2
        {
          ParticleCache* pc = new ParticleCacheDensitySelfContribution(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((ParticleCacheDensitySelfContribution*) pc)->stage() == m_stage);
          assert(((ParticleCacheDensitySelfContribution*) pc)->m_colour == m_colour);
          assert(((ParticleCacheDensitySelfContribution*) pc)->m_offset == m_offset);

          Particle::registerCache(pc);
          Particle::registerCache_0(this);
        }
      }
      // No? Then make a copy
      else 
      {
        ParticleCache* pc = new ParticleCacheDensitySelfContribution(*this);
        assert(pc->mySymbolName() == m_symbolName);
        assert(((ParticleCacheDensitySelfContribution*) pc)->stage() == m_stage);
        assert(((ParticleCacheDensitySelfContribution*) pc)->m_colour == m_colour);
        assert(((ParticleCacheDensitySelfContribution*) pc)->m_offset == m_offset);

        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);
          
          pc = new ParticleCacheDensitySelfContribution(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((ParticleCacheDensitySelfContribution*) pc)->stage() == m_stage);
          assert(((ParticleCacheDensitySelfContribution*) pc)->m_colour == m_colour);
          assert(((ParticleCacheDensitySelfContribution*) pc)->m_offset == m_offset);
          
          Particle::registerCache_0(pc);
        }
      
      }
    }
  }
  else /*it is a Symbol limited to one colour*/
  {
    m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
    
    // at least, if the attribute alrteady exists, we preserve the persistency and can hope that an other module has set it correctly
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
    {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype)
        throw gError("ParticleCacheDensitySelfContribution::setup", "Symbol " + m_symbolName + " already exists as a non-scalar.");
      else m_offset = Particle::s_tag_format[m_colour]./*indexOf*/offsetByName(m_symbolName);
    }
    else
      // if it does not yet exist we set the persistency to false
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;
    
    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else
    {
      ParticleCache* pc = new ParticleCacheDensitySelfContribution(*this);
      assert(pc->mySymbolName() == m_symbolName);
      assert(((ParticleCacheDensitySelfContribution*) pc)->stage() == m_stage);
      assert(((ParticleCacheDensitySelfContribution*) pc)->m_colour == m_colour);
      assert(((ParticleCacheDensitySelfContribution*) pc)->m_offset == m_offset);

      Particle::registerCache(pc);
      Particle::registerCache_0(this);
      
    }
  }
}

void ParticleCacheDensitySelfContribution::registerWithParticle()
{
}

