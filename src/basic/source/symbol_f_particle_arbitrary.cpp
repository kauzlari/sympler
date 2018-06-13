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

#include "symbol_f_particle_arbitrary.h"

#include <limits>

#define M_CONTROLLER ((Simulation*) m_parent) -> controller()


SymbolFParticleArbitrary::SymbolFParticleArbitrary(Simulation* parent)
  : ParticleCacheArbitrary(parent)
{
  init();

  // FIXME: preliminary (2018-05-24) hack for letting a unittest for
  // protected and hence not directly callable
  // SymbolFParticleArbitrary::setupOffset() pass. The only possible call via
  // SymbolFParticleArbitrary::setup() does not yet work in a unittest since it
  // throws an exception due to unset xml-input.
  m_offset = std::numeric_limits<std::size_t>::max();
}

SymbolFParticleArbitrary::~SymbolFParticleArbitrary()
{
}

void SymbolFParticleArbitrary::init()
{
  m_properties.setClassName("SymbolFParticleArbitrary");

  // Here we define the piece common to all children
  m_properties.setDescription
    (", which "
     "depends exclusively on other properties of the same particle. "
     "It may be defined by the attribute 'expression'. \nNote that, "
     "for this module, attributes 'overwrite' and 'symbol' must stay "
     "at their default values (see below). \nYou can only SET the "
     "force to a new value given by 'expression'. Incrementing the "
     "previous value is currently (2018-05-17) NOT possible. "
     "Therefore, two SymbolFParticleArbitrary will NOT obey the "
     "superposition principle for forces. You simply have to include "
     "everything in one module.");
  
  // It does not really matter, but since we have this option, and
  // since only 'true' makes sense when modifying the Particle
  // force, we overwrite the undefined setting from the parent
  // class by true
  m_overwrite = true;

  // The intention is to overwrite in the children the PREDEFINED setting from here
  // by a meaningful symbol name depending on the child class
  // SymbolFParticleArbitrary::checkISymbolNameUntouched() checks if user sets it to different
  // string and throws exception in that case. 
  m_symbolName = "PREDEFINED";

}

// experience shows that calling this method in
// SymbolFParticleArbitrary::setup() (which one could think would be
// ideal) is very unittest-unfriendly. So we call it in the children's
// setup() functions (if we do not forget, which is the drawback of
// this solution) 
void SymbolFParticleArbitrary::checkIfSymbolNameUntouched() const {

  if(m_symbolName != "PREDEFINED")
    throw gError
      ("SymbolFParticleArbitrary::checkIfSymbolNameUntouched: for module with className " + className() + " and name " + myName(), "Attribute 'symbol' was altered "
       "to " + m_symbolName + ". SymbolFParticleArbitrary does not permit "
       "to change this attribute.");
}

void SymbolFParticleArbitrary::setup()
{
  
  if(!m_overwrite)
    throw gError
      ("SymbolFParticleArbitrary::setup for module with className " + className() + " and name " + myName(), "Attribute 'overwrite' has "
       "forbidden value \"false\". SymbolFParticleArbitrary must have "
       "overwrite = \"true\"!");
  
  // following is for setting m_offset to the correct force_index in
  // each time step
  if(m_phase == 0) M_CONTROLLER->registerForPrecomputation_0(this);
  if(m_phase == 1) M_CONTROLLER->registerForPrecomputation(this);
  if(m_phase == 2) {
    M_CONTROLLER->registerForPrecomputation(this);
    M_CONTROLLER->registerForPrecomputation_0(this);
  }

  // All the stuff from the parent can be used if we have done what we
  // have done before in SymbolFParticleArbitrary::setup()
  ParticleCacheArbitrary::setup();

}

// This ParticleCache computes on Particle::v, hence no usage of
// Particle::s_tag_format[m_colour] methods such as addAttribute(..)
void SymbolFParticleArbitrary::setupOffset()
{
  // precompute() should set this to the correct force index, and a
  // unittest should check that it then has a reasonable value
  // here we set it to the largest possible size_t
  m_offset = std::numeric_limits<std::size_t>::max();
}

