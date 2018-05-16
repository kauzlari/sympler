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


#include "particle_vels.h"

const SymbolRegister<ParticleVels> particle_vels("ParticleVels");


ParticleVels::ParticleVels(Simulation* parent)
  : ParticleCacheArbitrary(parent)
{
  setFunctionReturnType();
  m_datatype = DataFormat::POINT;

  init();

  // FIXME: preliminary (2018-04-23) hack for letting a unittest for
  // protected and hence not directly callable
  // ParticleVels::setupOffset() pass. The only possible call via
  // ParticleVels::setup() dos not yet work in a unittest since it
  // throws and exception due to unset xml-input.
  m_offset = HUGE_VAL;
}

ParticleVels::~ParticleVels()
{
}

void ParticleVels::init()
{
  m_properties.setClassName("ParticleVels");
  m_properties.setDescription
    ("User defined modification of the particle velocity, which "
     "depends exclusively on other properties of the same particle. "
     "It may be defined by the attribute 'expression'. \nNote that, "
     "for this module, attributes 'overwrite' and 'symbol' must stay "
     "at their default values (see below). \nIf you want to add a "
     "term \"[newTerm]\" to the current velocity, then write "
     "'expression = \"[v]+[newTerm]\"'.");

  // It does not really matter, but since we have this option, and
  // since only 'true' makes sense when modifying the Particle
  // velocity, we overwrite the undefined setting from the parent
  // class by true
  m_overwrite = true;

  // The intention is to overwrite the undefined setting from the parent
  // class by the symbol name of the velocity.
  // ParticleVels::setup() checks if user sets it to different
  // string and throws exception in that case. 
  m_symbolName = "v";

}

void ParticleVels::setup()
{
  // Peculiarities of this class
  if(m_symbolName != "v")
    throw gError("ParticleVels::setup", "Attribute 'symbol' has forbidden value " + m_symbolName + ". ParticleVels must have symbol = \"v\"!");
  
  if(!m_overwrite)
    throw gError("ParticleVels::setup", "Attribute 'overwrite' has forbidden value \"false\". ParticleVels must have overwrite = \"true\"!");
  
  // All the stuff from the parent can be used if we have done what we
  // have done before in ParticleVels::setup()
  ParticleCacheArbitrary::setup();
}

// This ParticleCache computes on Particle::v, hence no usage of
// Particle::s_tag_format[m_colour] methods such as addAttribute(..)
void ParticleVels::setupOffset()
{
  // Should not be used. If used, the huge value hopefully produces at
  // least a seg.-fault
  m_offset = HUGE_VAL;

}
