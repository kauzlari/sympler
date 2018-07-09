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


#include "symbol_f_particle_scalar.h"

#include <limits>

#define M_CONTROLLER ((Simulation*) m_parent) -> controller()

const SymbolRegister<SymbolFParticleScalar> symbol_f_particle_scalar("SymbolFParticleScalar");


SymbolFParticleScalar::SymbolFParticleScalar(Simulation* parent)
  : SymbolFParticleArbitrary(parent)
{
  setFunctionReturnType();
  m_datatype = DataFormat::DOUBLE;

  init();
}

SymbolFParticleScalar::~SymbolFParticleScalar()
{
}

void SymbolFParticleScalar::init()
{
  m_properties.setClassName("SymbolFParticleScalar");

  // Add specific piece of description for this child class
  m_properties.setDescription
    (string("User defined modification of the newest flux on the "
	    "scalar given by attribute 'variable'")
     + m_properties.description()
     );

  STRINGPC
    (variable, m_variableName, "Name of the scalar variable for which "
     "the force/flux should be computed.");

  m_variableName = "UNDEFINED";  
}

void SymbolFParticleScalar::setup()
{

  // SymbolFParticleArbitrary implements this feature
  checkIfSymbolNameUntouched();
    
  // the most meaningful name we can assign in this class
  // This is done after
  // SymbolFParticleArbitrary::checkIfSymbolNameUntouched(), since the
  // latter still checks m_symbolName == "PREDEFINED", which is still
  // somewhat meaningful to do
  m_symbolName = string("force_" + m_variableName);
  
  // All the stuff from the parent can be used if we have done what we
  // have done before in SymbolFParticleScalar::setup()
  SymbolFParticleArbitrary::setup();

  // no reason why this couldn't be done at the very end. If run at the
  // top of this function, it prevents some unittests
  if(m_variableName == "UNDEFINED")
    throw gError
      ("SymbolFParticleScalar::setup", "Attribute 'variable' was not "
       "defined.");
}

// This ParticleCache computes on Particle::v, hence no usage of
// Particle::s_tag_format[m_colour] methods such as addAttribute(..)
void SymbolFParticleScalar::setupOffset()
{
  SymbolFParticleArbitrary::setupOffset();

  // with this we have each of the forces accessible and can assign it
  // to m_offset as needed (in precompute())
  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    // Currently (2018-05-24), setupOffset() is called by
    // ParticleCacheArbitrary::setup() in a context where m_colour has
    // the correct value.
    // FIXME: Is this an unsafe assumption? If yes, think about a safer
    // way.
    m_forceOffset[i] =
        Particle::s_tag_format[m_colour].attrByName(string("force_"
        + m_variableName + "_" + ObjToString(i))).offset;
  } 
}

