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


#include "particle_cache_arb_rng.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation *) m_parent)


ParticleCacheArbRNG::ParticleCacheArbRNG(Simulation* parent)
  : ParticleCacheArbitrary(parent)
{
  init();
}

ParticleCacheArbRNG::~ParticleCacheArbRNG()
{
}

void ParticleCacheArbRNG::init()
{

  INTPC(pLimit, m_plimit, INT_MIN, "Specifies the minimum number of particles for which the average of the random numbers of each time step will be shifted to zero while preserving the original variance.");

  m_plimit = 2;

}

void ParticleCacheArbRNG::setup()
{

  ParticleCacheArbitrary::setup();

  if(m_plimit < 1)
    throw gError("ParticleCacheArbRNG::setup", "For module " + m_properties.className() + ": Attribute 'pLimit' must have integer value > 0");

  m_phasePointer = M_SIMULATION->phase();

  // If this fails, we have to assign the phase in a setupAfterParticleCreation();
  assert(m_phasePointer);

  if (M_SIMULATION->randomize())
    m_rng.setSeed(getpid());
  else
    m_rng.setSeed(RNG_DEFAULT_SEED);

}

