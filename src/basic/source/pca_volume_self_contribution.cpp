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


#include "pca_volume_self_contribution.h"
#include "simulation.h"
#include "manager_cell.h"
#include "val_calculator_rho.h"
// #include "colour_pair.h"

const SymbolRegister<ParticleCacheVolumeSelfContribution> pca_volume_self("ParticleCacheVolumeSelfContribution");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

ParticleCacheVolumeSelfContribution::ParticleCacheVolumeSelfContribution
  (size_t colour, size_t offset, string symbolName, size_t density_offset, WeightingFunction *wf)
  : ParticleCache(colour, offset, symbolName)/*, m_offset(offset)*/, m_density_offset(density_offset), m_wf(wf)
{
  // depends on loacal density (m_stage = 0)
  m_stage = 1;
  m_datatype = DataFormat::DOUBLE;
}

ParticleCacheVolumeSelfContribution::ParticleCacheVolumeSelfContribution(/*Node*/Simulation* parent)
  : ParticleCache(parent)
{
  m_stage = 1;
  m_datatype = DataFormat::DOUBLE;
  init(); 
}

ParticleCacheVolumeSelfContribution::~ParticleCacheVolumeSelfContribution()
{
}

void ParticleCacheVolumeSelfContribution::init()
{
  m_properties.setClassName("ParticleCacheVolumeSelfContribution");
  m_properties.setName("VolumeSC");

  m_properties.setDescription("Contribution of the particle itself to its local volume.");
      
  STRINGPC
      (symbol, m_symbolName,
       "Name for the local volume.");
  
  STRINGPC
      (density, m_densitySymbol,
       "Name of the symbol to take the local density from.");
  
  STRINGPC
      (weightingFunction, m_wfName,
       "Defines the weighting function to be used for the computation of the local density.");
     
  m_wfName = "default";
  m_symbolName = "V";
  m_densitySymbol = "undefined";   
}

void ParticleCacheVolumeSelfContribution::setup()
{
  if(m_species == "undefined")
    throw gError("ParticleCacheVolumeSelfContribution::setup", "Attribute 'species' has value \"undefined\""); 

  if(m_densitySymbol == "undefined")
    throw gError("ParticleCacheVolumeSelfContribution::setup", "Attribute 'volumeSymbol' has value \"undefined\""); 

  m_wf = M_SIMULATION->findWeightingFunction(m_wfName);

  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours too?
  if(m_species == "ALL")
  {
    // is the symbol already existing somewhere?
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
        throw gError("ParticleCacheVolumeSelfContribution::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed with this Symbol calculator");
    }

    bool newCache = true;
    m_colour = 0;
      
    FOR_EACH_COLOUR_PAIR
    (
      M_MANAGER,
      
      // ValCalculatorRho decides about persistency for density symbol
    if(m_phaseUser == 0)
      cp->registerCalc_0(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), true/*oneProp*/);
    if(m_phaseUser == 1)
      cp->registerCalc(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), true/*oneProp*/);
    if(m_phaseUser == 2)
    {
      cp->registerCalc_0(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), true/*oneProp*/);
      cp->registerCalc(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), true/*oneProp*/);
    }
//       bool newCache = false;
      
      if(m_colour < cp->firstColour() || m_colour < cp->secondColour())
      {
        ++m_colour;
        newCache = true;
      }
      // if density is persistent, so is the volume
      bool persistency = Particle::s_tag_format[m_colour].attrByName(m_densitySymbol).persistent;
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, persistency, m_symbolName).offset;
      
      if(newCache)
      {
        if(cp->firstColour() == m_colour) m_density_offset = tempPair.first;
        else
        {
          assert(cp->secondColour() == m_colour);
          m_density_offset = tempPair.second;
        }
          // is it the last cache to be created?
        if(m_colour == M_MANAGER->nColours()-1)
        {
          if(m_phaseUser == 0)
            Particle::registerCache_0(this);
          else if(m_phaseUser == 1)
            Particle::registerCache(this);
          else // so it is 2
          {
            ParticleCache* pc = new ParticleCacheVolumeSelfContribution(*this);
            assert(pc->mySymbolName() == m_symbolName);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_densitySymbol == m_densitySymbol);
            assert(pc->stage() == m_stage);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_colour == m_colour);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_offset == m_offset);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_density_offset == m_density_offset);
  
            Particle::registerCache(pc);
            Particle::registerCache_0(this);
          }
        }
          // No? Then make a copy
        else 
        {
          ParticleCache* pc = new ParticleCacheVolumeSelfContribution(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((ParticleCacheVolumeSelfContribution*) pc)->m_densitySymbol == m_densitySymbol);
          assert(pc->stage() == m_stage);
          assert(((ParticleCacheVolumeSelfContribution*) pc)->m_colour == m_colour);
          assert(((ParticleCacheVolumeSelfContribution*) pc)->m_offset == m_offset);
          assert(((ParticleCacheVolumeSelfContribution*) pc)->m_density_offset == m_density_offset);
  
          if(m_phaseUser == 0)
            Particle::registerCache_0(pc);
          else if(m_phaseUser == 1)
            Particle::registerCache(pc);
          else // so it is 2
          {
            Particle::registerCache_0(pc);
            
            pc = new ParticleCacheVolumeSelfContribution(*this);
            assert(pc->mySymbolName() == m_symbolName);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_densitySymbol == m_densitySymbol);
            assert(pc->stage() == m_stage);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_colour == m_colour);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_offset == m_offset);
            assert(((ParticleCacheVolumeSelfContribution*) pc)->m_density_offset == m_density_offset);
             
            Particle::registerCache(pc);
          }
        }
      }
    );

  }
  else /*it is a Symbol limited to one colour*/
  {
    m_colour = M_MANAGER->getColour/*AndAdd*/(m_species); 
    // is the Symbol already existing?
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      throw gError("ParticleCacheDensity0Oc::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed with this Symbol calculator");
    
    if(m_densitySymbol == "undefined")
      throw gError("ParticleCacheVolumeSelfContribution::setup", "Attribute 'volumeSymbol' has value \"undefined\""); 

    // only the homogeneous ColourPair will contribute
    ColourPair* cp = M_MANAGER->cp(m_colour, m_colour);
    
    // ValCalculatorRho decides about persistency for density symbol
    if(m_phaseUser == 0)
      cp->registerCalc_0(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), false/*oneProp*/);
    if(m_phaseUser == 1)
      cp->registerCalc(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), false/*oneProp*/);
    if(m_phaseUser == 2)
    {
      cp->registerCalc_0(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), false/*oneProp*/);
      cp->registerCalc(tempPair, new ValCalculatorRho(/*m_wf, */m_densitySymbol, m_wf/*, extra_string*//*, false*//*persistency*/), false/*oneProp*/);
    }
    m_density_offset = tempPair.first;
    
    // if density symbol is persistent, so is the volume symbol
    bool persistency = Particle::s_tag_format[m_colour].attrByName(m_densitySymbol).persistent;
    m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, persistency, m_symbolName).offset;
    
    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else // so it is 2
    {
      ParticleCache* pc = new ParticleCacheVolumeSelfContribution(*this);
      assert(pc->mySymbolName() == m_symbolName);
      assert(((ParticleCacheVolumeSelfContribution*) pc)->m_densitySymbol == m_densitySymbol);
      assert(pc->stage() == m_stage);
      assert(((ParticleCacheVolumeSelfContribution*) pc)->m_colour == m_colour);
      assert(((ParticleCacheVolumeSelfContribution*) pc)->m_offset == m_offset);
      assert(((ParticleCacheVolumeSelfContribution*) pc)->m_density_offset == m_density_offset);

      Particle::registerCache_0(pc);
      Particle::registerCache(this);
            
    }
  }

}

void ParticleCacheVolumeSelfContribution::registerWithParticle()
{
}

