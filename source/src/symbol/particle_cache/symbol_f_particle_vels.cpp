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
  : SymbolFParticleArbitrary(parent)
{
  setFunctionReturnType();
  m_datatype = DataFormat::POINT;

  init();
}

SymbolFParticleVels::~SymbolFParticleVels()
{
}

void SymbolFParticleVels::init()
{
  m_properties.setClassName("SymbolFParticleVels");
  
  m_properties.setDescription
    (string("User defined modification of the newest particle force")
     + m_properties.description()
     );
}

void SymbolFParticleVels::setup()
{
  // SymbolFParticleArbitrary implements this feature
  checkIfSymbolNameUntouched();
  
  // The most meaningful we can assign (?)
  m_symbolName = "__force__";
    
  // All the stuff from the parent can be used if we have done what we
  // have done before in SymbolFParticleVels::setup()
  SymbolFParticleArbitrary::setup();
}

// This ParticleCache computes on Particle::v, hence no usage of
// Particle::s_tag_format[m_colour] methods such as addAttribute(..)
void SymbolFParticleVels::setupOffset()
{
  SymbolFParticleArbitrary::setupOffset();  
}

