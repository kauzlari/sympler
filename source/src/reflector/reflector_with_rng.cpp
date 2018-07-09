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



#include "reflector_with_rng.h"
#include "manager_cell.h"
#include "simulation.h"

// #include "valgrind/memcheck.h"

#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation *) M_PHASE->parent())

#define M_CONTROLLER M_SIMULATION->controller()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

ReflectorWithRng::ReflectorWithRng(Boundary* boundary): Reflector(boundary /*wall*/)
{
  init();
}


ReflectorWithRng::~ReflectorWithRng()
{
}


//---- Methods ----

void ReflectorWithRng::init()
{
  m_properties.setClassName("ReflectorWithRng");

  m_properties.setDescription(
    "Reflector with random number generator."
  );
}

void ReflectorWithRng::setup()
{
  Reflector::setup();

  if (M_SIMULATION->randomize()) {
    MSG_DEBUG("ReflectorWithRng::setup", "randomizing");
    m_rng.setSeed(getpid());
  }
  else {
    MSG_DEBUG("ReflectorWithRng::setup", "NOT randomizing and using seed " << RNG_DEFAULT_SEED);
    m_rng.setSeed(RNG_DEFAULT_SEED);
  }
}


