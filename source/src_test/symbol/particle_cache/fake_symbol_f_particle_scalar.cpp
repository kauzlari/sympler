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


#include "fake_symbol_f_particle_scalar.h"

FakeSymbolFParticleScalar::FakeSymbolFParticleScalar(Simulation* parent)
  : SymbolFParticleScalar(parent)
{
}

FakeSymbolFParticleScalar::~FakeSymbolFParticleScalar()
{
}

void FakeSymbolFParticleScalar::setupForceOffsetForTestsOnly()
{
  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i) {
    m_forceOffset[i] =
      Particle::s_tag_format[m_colour]
      .attrByName(string("TESTforce_" + m_variableName + "_"
			 + ObjToString(i))).offset;
  } 
}

void FakeSymbolFParticleScalar::setColour(size_t colour)
{
  m_colour = colour;
}
