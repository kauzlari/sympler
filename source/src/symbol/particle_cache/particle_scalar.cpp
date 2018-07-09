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


#include "particle_scalar.h"

const SymbolRegister<ParticleScalar> particle_scalar("ParticleScalar");


ParticleScalar::ParticleScalar(Simulation* parent)
  : ParticleCacheArbitrary(parent)
{
  setFunctionReturnType();
  m_datatype = DataFormat::DOUBLE;
  init();
}

ParticleScalar::~ParticleScalar()
{
}

void ParticleScalar::init()
{
  m_properties.setClassName("ParticleScalar");
  m_properties.setDescription("User defined scalar particle property, which depends exclusively on other properties of the same particle. It may be defined by the attribute 'expression'.");
}

void ParticleScalar::setup()
{
  ParticleCacheArbitrary::setup();
}

