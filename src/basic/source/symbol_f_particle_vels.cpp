/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#include "symbol_f_particle_vels.h"

#include <limits>

#define M_CONTROLLER ((Simulation*) m_parent) -> controller()

const SymbolRegister<SymbolFParticleVels> symbol_f_particle_vels("SymbolFParticleVels");


SymbolFParticleVels::SymbolFParticleVels(Simulation* parent)
  : ParticleCacheArbitrary(parent)
{
  setFunctionReturnType();
  m_datatype = DataFormat::POINT;

  init();

  // FIXME: preliminary (2018-05-17) hack for letting a unittest for
  // protected and hence not directly callable
  // SymbolFParticleVels::setupOffset() pass. The only possible call via
  // SymbolFParticleVels::setup() does not yet work in a unittest since it
  // throws an exception due to unset xml-input.
  m_offset = std::numeric_limits<std::size_t>::max();
}

SymbolFParticleVels::~SymbolFParticleVels()
{
}

void SymbolFParticleVels::init()
{
  m_properties.setClassName("SymbolFParticleVels");
  
  m_properties.setDescription
    ("User defined modification of the newest particle force, which "
     "depends exclusively on other properties of the same particle. "
     "It may be defined by the attribute 'expression'. \nNote that, "
     "for this module, attributes 'overwrite' and 'symbol' must stay "
     "at their default values (see below). \nYou can only SET the "
     "force to a new value given by 'expression'. Incrementing the "
     "previous value is currently (2018-05-17) NOT possible. "
     "Therefore, two SymbolFParticleVels will NOT obey the "
     "superposition principle for forces. You simply have to include "
     "everything in one module.");

  // It does not really matter, but since we have this option, and
  // since only 'true' makes sense when modifying the Particle
  // force, we overwrite the undefined setting from the parent
  // class by true
  m_overwrite = true;

  // The intention is to overwrite the undefined setting from the parent
  // class by a symbol name that is (hopefully) used nowhere else
  // SymbolFParticleVels::setup() checks if user sets it to different
  // string and throws exception in that case. 
  m_symbolName = "__force__";

}

void SymbolFParticleVels::setup()
{
  // Peculiarities of this class

  if(m_symbolName != "__force__")
    throw gError("SymbolFParticleVels::setup", "Attribute 'symbol' has forbidden value " + m_symbolName + ". SymbolFParticleVels must have default symbol name \"__force__\"!");
  
  if(!m_overwrite)
    throw gError("SymbolFParticleVels::setup", "Attribute 'overwrite' has forbidden value \"false\". SymbolFParticleVels must have overwrite = \"true\"!");

  // following is for setting m_offset to the correct force_index in
  // each time step
  if(m_phase == 0) M_CONTROLLER->registerForPrecomputation_0(this);
  if(m_phase == 1) M_CONTROLLER->registerForPrecomputation(this);
  if(m_phase == 2) {
    M_CONTROLLER->registerForPrecomputation(this);
    M_CONTROLLER->registerForPrecomputation_0(this);
  }
  
  // All the stuff from the parent can be used if we have done what we
  // have done before in SymbolFParticleVels::setup()
  ParticleCacheArbitrary::setup();
}

// This ParticleCache computes on Particle::v, hence no usage of
// Particle::s_tag_format[m_colour] methods such as addAttribute(..)
void SymbolFParticleVels::setupOffset()
{
  // precompute() should set this to the correct force index, and a
  // unittest should check that it then has a reasonable value
  // here we set it to the largest possible size_t
  m_offset = std::numeric_limits<std::size_t>::max();
  
}

