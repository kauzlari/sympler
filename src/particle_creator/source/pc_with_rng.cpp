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


/* $*/


#include <fstream>
#include <iostream>

#include "phase.h"
#include "manager_cell.h"
#include "simulation.h"


#include "pc_with_rng.h"

// #define M_BOUNDARY ((Boundary*) m_parent)
// #define M_PHASE ((Phase*) M_BOUNDARY->parent())
// #define M_SIMULATION ((Simulation*) M_PHASE->parent())
// #define M_MANAGER M_PHASE->manager()


/* ParticleCreatorWithRngPCalc */

ParticleCreatorWithRngPCalc::ParticleCreatorWithRngPCalc(Boundary *boundary) : ParticleCreatorFreePCalc(boundary)
{
  init();
}

ParticleCreatorWithRngPCalc::~ParticleCreatorWithRngPCalc()
{
}

//--- Methods ---

void ParticleCreatorWithRngPCalc::init()
{
  m_properties.setClassName("ParticleCreatorWithRngPCalc");

  /* Allow unknown properties. Those ones have to be identified later.
  They are used to set the particles degree of freedoms initially. */
  m_properties.allowUnknown();

  BOOLPC
      (randomize, m_randomize,
       "If set to \"yes\", this produces a different set of random numbers for the positions and velocities of the particles.");

  DOUBLEPC
    (kBToverM, m_temperature, 0,
     "k_BT/m = <v^2> for setting the initial thermal velocities of the particles. Note that this "
     "can be overriden by setting vel*.");

  m_randomize = false;
  m_temperature = 1;

}

void ParticleCreatorWithRngPCalc::setup()
{
  ParticleCreatorFreePCalc::setup();


  if (m_randomize)
    m_rng.setSeed(getpid());
  else
    m_rng.setSeed(RNG_DEFAULT_SEED);


}

/* ParticleCreatorWithRngF */

ParticleCreatorWithRngF::ParticleCreatorWithRngF(Boundary *boundary) : ParticleCreatorFreeF(boundary)
{
  init();
}

ParticleCreatorWithRngF::~ParticleCreatorWithRngF()
{
}

//--- Methods ---

void ParticleCreatorWithRngF::init()
{
  m_properties.setClassName("ParticleCreatorWithRngF");

  /* Allow unknown properties. Those ones have to be identified later.
  They are used to set the particles degree of freedoms initially. */
  m_properties.allowUnknown();

  BOOLPC
      (randomize, m_randomize,
       "If set to \"yes\", this produces a different set of random numbers for the positions and velocities of the particles.");

  m_randomize = false;
}

void ParticleCreatorWithRngF::setup()
{
  ParticleCreatorFreeF::setup();

  if (m_randomize)
    m_rng.setSeed(getpid());
  else
    m_rng.setSeed(RNG_DEFAULT_SEED);

}

