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


#include "pca_density_0oc.h"
#include "colour_pair.h"
#include "val_calculator_volume.h"

const SymbolRegister<ParticleCacheDensity0Oc> pca_density_0oc(/*"ParticleCacheDensity0Oc"*/"Density0Oc");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

ParticleCacheDensity0Oc::ParticleCacheDensity0Oc
  (/*Node*/Simulation* parent/*ColourPair* cp_arg, size_t colour, size_t density_offset, WeightingFunction *wf, string extra_string,  bool oneProp*/)
  : ParticleCache(parent)/*ParticleCache(colour), m_density_offset(NULL), m_wf(wf)*/
{
  // depends on local volume (m_stage = 1), which again depends on local 
  // density (m_stage = 0)
  m_stage = 2;
  m_datatype = DataFormat::DOUBLE;
  init();
    
  // next two lines should be redundant, since also in ValCalculatorVolume::SetSlots,
  // but they were not commented out until the implementation of Symbol hierarchie
/*  cp_arg->setCutoff(wf->cutoff());
  cp_arg->setNeedPairs(true);*/

  
/*  pair<size_t, size_t> tempPair;
  // create a ValCalculatorVolume
  // this should place the PCAVolume in front of the PCADensity0Oc
        // first we register the value in the CPs != m_cp, if m_oneProp = true
  
  MSG_DEBUG("ParticleCacheDensity0Oc::ParticleCacheDensity0Oc", "creating ValCalculatorVolume now");
  if(oneProp)
  {
    FOR_EACH_COLOUR_PAIR
        (
        cp_arg->manager(),
        // if this is m_cp then do nothing (will be done afterwards)
      if(cp_arg->firstColour() != cp->firstColour() || cp_arg->secondColour() != cp->secondColour())
      {
        cp->registerCalc(tempPair, new ValCalculatorVolume(m_wf, extra_string), true);
      }
        );
  }
  cp_arg->registerCalc(tempPair, new ValCalculatorVolume(wf, extra_string), oneProp);
  if(cp_arg->firstColour() == colour) m_volume_offset = tempPair.first;
  else
  {
    assert(cp_arg->secondColour() == colour);
    m_volume_offset = tempPair.second;
  }*/
}


ParticleCacheDensity0Oc::~ParticleCacheDensity0Oc()
{
}

void ParticleCacheDensity0Oc::registerWithParticle()
{
}

void ParticleCacheDensity0Oc::init()
{
  m_properties.setClassName("Density0Oc");

  m_properties.setDescription("0th order correction for the local density \"n_c = n/V\". \"n_c\" is the corrected density. \"n\" is the local density computed through \"n_i\" = Sum(j, m_j*W_ij) for a particle i and all its neighbours j. \"m\" is the particle mass and W_ij is a weighting function. \"V\" is the local volume computed through \"V_i = Sum(j, W_ij/n_j)\" for a particle i and all its neighbours j. This is used, e.g., in J. Non-Newt. Fluid Mech. 132, p.61.\nThe Symbol \"V\" is additionally introduced by this Calculator.");
      
  STRINGPC
      (symbol, m_symbolName,
       "Name for the local density.");
  
  STRINGPC(volumeSymbol, m_volumeSymbolName, "Name to be given to the Symbol for the local volume of a particle");

  STRINGPC
      (weightingFunction, m_weighting_function,
       "Defines the weighting function to be used for the computation of the local density and the local volume.");
     
  m_weighting_function = "default";
  
  m_volumeSymbolName = "undefined";
}

void ParticleCacheDensity0Oc::setup()
{
  if(m_species == "undefined")
    throw gError("ParticleCacheDensity0Oc::setup", "Attribute 'species' has value \"undefined\""); 
  
  if(m_volumeSymbolName == "undefined")
    throw gError("ParticleCacheDensity0Oc::setup", "Attribute 'volumeSymbol' has value \"undefined\""); 
  
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("ParticleCacheDensity0Oc::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  
  m_wf = M_SIMULATION->findWeightingFunction(m_weighting_function);
  
  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours too?
  if(m_species == "ALL")
  {
    // is the symbol already existing somewhere?
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
        throw gError("ParticleCacheDensity0Oc::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed with this Symbol calculator");
    }

    m_colour = 0;
    bool newCache = true;  
    
    // FIXME: is there a check for the kind of persistency missing? Currently it is always false. Must it be true if there is no IntegratorPosition?
//    m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName);
    
    FOR_EACH_COLOUR_PAIR
    (
      M_MANAGER,
//       cp->registerCalc(tempPair, new ValCalculatorRho(m_wf, extra_string), true);
      
      // ValCalculatorVolume should create ValCalculatorRho and those two should take care for persistency
    if(m_phaseUser == 0)
      cp->registerCalc_0(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), true/*oneProp*/);
    if(m_phaseUser == 1)
      cp->registerCalc(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), true/*oneProp*/);
    if(m_phaseUser == 2)
    {
      cp->registerCalc_0(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), true/*oneProp*/);
      cp->registerCalc(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), true/*oneProp*/);
    } 
//       bool newCache = false;
      
      if(m_colour < cp->firstColour() || m_colour < cp->secondColour())
      {
        ++m_colour;
        newCache = true;
      }

      // since the ValCalcs are already created, the attribute should already exist
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype)
        throw gError("ParticleCacheDensity0Oc::setup", "Symbol '" + m_symbolName + "' already exists as a non-scalar.");
      m_offset = Particle::s_tag_format[m_colour].offsetByName/*indexOf*/(m_symbolName/*, m_datatype, false*//*persistency*//*, m_symbolName*/);
    
      if(newCache)
      {
        if(cp->firstColour() == m_colour) m_volume_offset = tempPair.first;
        else
        {
          assert(cp->secondColour() == m_colour);
          m_volume_offset = tempPair.second;
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
            ParticleCache* pc = new ParticleCacheDensity0Oc(*this);
            assert(pc->mySymbolName() == m_symbolName);
            assert(((ParticleCacheDensity0Oc*) pc)->myVolumeSymbolName() == m_volumeSymbolName);
            assert(pc->stage() == m_stage);
            assert(((ParticleCacheDensity0Oc*) pc)->m_colour == m_colour);
            assert(((ParticleCacheDensity0Oc*) pc)->m_offset == m_offset);
            assert(((ParticleCacheDensity0Oc*) pc)->m_volume_offset == m_volume_offset);

            Particle::registerCache(pc);
            
            Particle::registerCache_0(this);
          }
        }
        // No? Then make a copy
        else 
        {
          ParticleCache* pc = new ParticleCacheDensity0Oc(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((ParticleCacheDensity0Oc*) pc)->myVolumeSymbolName() == m_volumeSymbolName);
          assert(pc->stage() == m_stage);
          assert(((ParticleCacheDensity0Oc*) pc)->m_colour == m_colour);
          assert(((ParticleCacheDensity0Oc*) pc)->m_offset == m_offset);
          assert(((ParticleCacheDensity0Oc*) pc)->m_volume_offset == m_volume_offset);

          if(m_phaseUser == 0)
            Particle::registerCache_0(pc);
          else if(m_phaseUser == 1)
            Particle::registerCache(pc);
          else // so it is 2
          {
            Particle::registerCache_0(pc);
            
            pc = new ParticleCacheDensity0Oc(*this);
            assert(pc->mySymbolName() == m_symbolName);
            assert(((ParticleCacheDensity0Oc*) pc)->myVolumeSymbolName() == m_volumeSymbolName);
            assert(pc->stage() == m_stage);
            assert(((ParticleCacheDensity0Oc*) pc)->m_colour == m_colour);
            assert(((ParticleCacheDensity0Oc*) pc)->m_offset == m_offset);
            assert(((ParticleCacheDensity0Oc*) pc)->m_volume_offset == m_volume_offset);
            
            Particle::registerCache(pc);
          }
          
          newCache = false;
      
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
    
    if(m_volumeSymbolName == "undefined")
      throw gError("ParticleCacheDensity0Oc::setup", "Attribute 'volumeSymbol' has value \"undefined\""); 
  
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
    {
      // with the following, we preserve the existing persistency status
      if(Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype == DataFormat::DOUBLE)
        m_offset = Particle::s_tag_format[m_colour].offsetByName/*indexOf*/(m_symbolName);
      else 
        throw gError("ParticleCacheDensity0Oc::setup", "Symbol + " + m_symbolName + " already exists as a non-scalar.");
    }
    else
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;
    
    // only the homogeneous ColourPair will contribute
    ColourPair* cp = M_MANAGER->cp(m_colour, m_colour);
    
    
    if(Particle::s_tag_format[m_colour].attrExists(m_volumeSymbolName))
    {
      // with the following, we preserve the existing persistency status
      if(Particle::s_tag_format[m_colour].attrByName(m_volumeSymbolName).datatype == DataFormat::DOUBLE)
        m_volume_offset = Particle::s_tag_format[m_colour].offsetByName/*indexOf*/(m_volumeSymbolName);
      else throw gError("ParticleCacheDensity0Oc::setup", "Symbol + " + m_volumeSymbolName + " already exists as a non-scalar.");
    }
    else
    {
      
      // ValCalculatorVolume should create ValCalculatorRho and take care for persistency
      if(m_phaseUser == 0)
        cp->registerCalc_0(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), false/*oneProp*/);
      if(m_phaseUser == 1)
        cp->registerCalc(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), false/*oneProp*/);
      m_volume_offset = tempPair.first;
      if(m_phaseUser == 2)
      {
        cp->registerCalc_0(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), false/*oneProp*/);
        cp->registerCalc(tempPair, new ValCalculatorVolume(/*m_wf, */m_symbolName, m_volumeSymbolName/*, extra_string*//*, false*//*persistency*/, m_wf), false/*oneProp*/);
      }
    }
    
    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else // so it is 2
    {
      ParticleCache* pc = new ParticleCacheDensity0Oc(*this);
      assert(pc->mySymbolName() == m_symbolName);
      assert(((ParticleCacheDensity0Oc*) pc)->myVolumeSymbolName() == m_volumeSymbolName);
      assert(pc->stage() == m_stage);
      assert(((ParticleCacheDensity0Oc*) pc)->m_colour == m_colour);
      assert(((ParticleCacheDensity0Oc*) pc)->m_offset == m_offset);
      assert(((ParticleCacheDensity0Oc*) pc)->m_volume_offset == m_volume_offset);
            
      Particle::registerCache_0(pc);
      Particle::registerCache(this);
      
    }
  }

}
