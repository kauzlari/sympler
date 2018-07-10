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



#include "pca_energy_entropy.h"
// #include "simulation.h"
// #include "manager_cell.h"

// const SymbolRegister<ParticleCacheEnergyEntropy> pca_energy_entropy("ParticleCacheEnergyEntropy");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

// the symbolName is useless for this Cache
// the offset is set in registerWithParticle
ParticleCacheEnergyEntropy::ParticleCacheEnergyEntropy
  (size_t colour/*, size_t offset, string symbolName*/, double density_cutoff, ColourPair *cp, WeightingFunction *wf,
   FunctionFixed *de_dn, FunctionFixed *Tds_dn)
  : ParticleCache(colour, 1000000/*offset*/, ""/*symbolName*/), m_de_dn(de_dn), m_Tds_dn(Tds_dn), m_cp(cp),
 m_wf(wf), m_density_cutoff(density_cutoff)
{
  // depends on loacal density (m_stage = 0)
//   m_colour = colour;
  m_stage = 1;
  m_datatype = DataFormat::DOUBLE;
//   m_symbolName = "";
  m_Tds_dn_symbol = "";
}

ParticleCacheEnergyEntropy::ParticleCacheEnergyEntropy(/*Node*/Simulation* parent)
  : ParticleCache(parent)
{
  m_stage = 1;
  m_datatype = DataFormat::DOUBLE;
//   init();
}

ParticleCacheEnergyEntropy::~ParticleCacheEnergyEntropy()
{
}

// the following does not work (but is not called by any constructor) 
// because the pointers to the functions are undefined
void ParticleCacheEnergyEntropy::init()
{
  m_properties.setClassName("ParticleCacheEnergyEntropy");
  
  m_properties.setDescription("Used to calculate the per particle energy and entropy.");
      
  FUNCTIONFIXEDPC(dedn, *m_de_dn, "Algebraic expression for de/dn\n, i.e., the derivative of the local configurational energy 'e' with respect to the local density 'n'. 'e' and 'n' are allowed variables.\nExample: For the van-der-Waals gas de/dn = -a*N, where N is the number of microscopic entities, the fluid particle consists of and a is the van-der-Waals parameter.");
  
  STRINGPC(dednSymbol, m_symbolName, "Symbol name of de/dn.");
  
  FUNCTIONFIXEDPC
      (Tds_dn, *m_Tds_dn,
       "Algebraic expression for T*ds/dn = -Pressure, i.e., the derivative of the local entropy with respect to the local density (times the temperature).  'e' and 'n' are allowed variables.\nExample: For the van-der-Waals gas T*ds/dn = -kB*T*N/(n*(1-n*b)), if we ignore long range contributions in the fluid particles' entropy. T is temperature, s is local entropy, n is local density, kB is the Boltzmann constant, N is the number of microscopic entities, the fluid particle consists of and b is the van-der-Waals parameter.");

  STRINGPC(TdsdnSymbol, m_Tds_dn_symbol, "Symbol name of the derivative of the local entropy with respect to the local density (times the temperature).");

  m_de_dn->addVariables("n", "e");  
  m_Tds_dn->addVariables("n", "e");
    
  m_symbolName = "undefined";
  m_Tds_dn_symbol = "undefined";
}

void ParticleCacheEnergyEntropy::setup()
{
  if(m_species == "undefined")
    throw gError("ParticleCacheEnergyEntropy::setup", "Attribute 'species' has value \"undefined\""); 
  if(m_species == "ALL")
    throw gError("ParticleCacheEnergyEntropy::setup", "Sorry, this Calculator cannot yet handle the option \"ALL\" for 'species'."); 

  if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
    throw gError("ParticleCacheEnergyEntropy::setup", "Symbol " + m_symbolName + " is already existing for species '" + m_species + "'. Second definition is not allowed with this Symbol calculator");
  if(Particle::s_tag_format[m_colour].attrExists(m_Tds_dn_symbol))
    throw gError("ParticleCacheEnergyEntropy::setup", "Symbol " + m_Tds_dn_symbol + " is already existing for species '" + m_species + "'. Second definition is not allowed with this Symbol calculator");
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("ParticleCacheEnergyEntropy::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");

  m_cp = M_MANAGER->cp(m_colour, m_colour);
  
  // next will also call ParticleCacheEnergyEntropy::registerWithParticle
  if(m_phaseUser == 0)
    Particle::registerCache_0(this);
  if(m_phaseUser == 1)
    Particle::registerCache(this);
  else // so it is 2
  {
    throw gError("ParticleCacheEnergyEntropy::setup", "Sorry, this Calculator cannot yet handle the option \"2\" for 'stage'."); 
  }
}

void ParticleCacheEnergyEntropy::registerWithParticle()
{
  pair<size_t, size_t> dummy;

  assert(m_cp->firstColour() == m_cp->secondColour());

  // the following is commented out because I think it does not 
  // make sense for a Cache that is limited to one colour

/*        // first we register the value in the CPs != m_cp, if m_oneProp = true
  FOR_EACH_COLOUR_PAIR
      (
      m_cp->manager(),
      // if this is m_cp then do nothing (will be done afterwards)
  if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
  {
    cp->registerCalc(dummy, new ValCalculatorRho(m_wf, ""), true);
  }
      );*/
  
  if(m_phaseUser == 0)
    m_cp->registerCalc_0(dummy, new ValCalculatorRho(m_densitySymbol, m_wf/*, ""*/), true);
  if(m_phaseUser == 1)
    m_cp->registerCalc(dummy, new ValCalculatorRho(m_densitySymbol, m_wf/*, ""*/), true);
  if(m_phaseUser == 2)
    throw gError("ParticleCacheEnergyEntropy::setup", "Sorry, this Calculator cannot yet handle the option \"2\" for 'stage'."); 
      
  /*m_density_o*/m_offset = dummy.first;

  m_energy_o = Particle::tagFormat(m_colour).attrByName("internal_energy").offset;

  if(m_symbolName == "") m_symbolName = string("de_dn_") + m_cp->toString() + "_" + m_wf->name();
  
  m_de_dn_o = Particle::tagFormat(m_colour).addAttribute
      (m_symbolName, DataFormat::DOUBLE).offset;
  
  if(m_Tds_dn_symbol == "") m_Tds_dn_symbol = string("Tds_dn_") + m_cp->toString() + "_" + m_wf->name();

  m_Tds_dn_o = Particle::tagFormat(m_colour).addAttribute
      (m_Tds_dn_symbol, DataFormat::DOUBLE).offset;

  assert(m_offset < 1000000);
}

