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


#include "particle_rand_norm_vector.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER (M_SIMULATION->controller())

const SymbolRegister<ParticleRandNormVector> particle_rand_norm_vector("ParticleRandNormVector");


ParticleRandNormVector::ParticleRandNormVector(Simulation* parent)
  : ParticleCacheArbRNG(parent)
{
  setFunctionReturnType();
  m_datatype = DataFormat::POINT;
  init();
}

ParticleRandNormVector::~ParticleRandNormVector()
{
}

void ParticleRandNormVector::init()
{
  m_properties.setClassName("ParticleRandNormVector");
  m_properties.setDescription("User-defined symbol for a particle computing a vector of 3 statisticaly independent random normaly distributed scalars with unit variance multiplied by previously computed particle properties (symbols). The latter may be defined by the attribute 'expression'. The 'expression' can also be used to modify the variance of the random number (by including the square root of the desired variance). This module's seed is randomised by the attribute 'randomize' of module 'Simulation'. The average of the ramdom numbers will be shifted to exactly zero for each individual time step, unless this is explicitly switched off by setting the attribute 'pLimit' larger than the number of particles. In the latter case, the average over an infinite number of time steps is still zero.");

  m_expression = "idVec(1)";

}

void ParticleRandNormVector::setup()
{
  ParticleCacheArbRNG::setup();

  if(m_phase == 0) M_CONTROLLER->registerForPrecomputation_0(this);
  if(m_phase == 1) M_CONTROLLER->registerForPrecomputation(this);
  if(m_phase == 2) {
    M_CONTROLLER->registerForPrecomputation(this);
    M_CONTROLLER->registerForPrecomputation_0(this);
  }
}

