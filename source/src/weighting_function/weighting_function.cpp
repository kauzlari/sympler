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



#include "weighting_function.h"

#include "simulation.h"

REGISTER_SMART_ENUM
(WeightingFunction_Factory,
 "Weighting functions define SPH kernels and their respective derivatives.");


#define M_SIMULATION ((Simulation*) m_parent)


/* --- WeightingFunction --- */

/* Constructor/Destructor */

WeightingFunction::WeightingFunction()
{
  throw gError("WeightingFunction::WeightingFunction(default)", "Should not be called. Contact the programmer.");
}


WeightingFunction::WeightingFunction(Node *parent): Node(parent) {
  init();
}

/*
WeightingFunction(double rc): m_rc(rc) {
  prepare();
}
*/

WeightingFunction::~WeightingFunction() {
}



/* Methods */

void WeightingFunction::init()
{
  m_properties.setClassName("WeightingFunction");

  STRINGPC
    (name, m_name,
     "Name for the weighting function.");

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Cut-off radius for the weighting function.");

  m_name = "default";
  m_cutoff = 1;
}



/* --- WeightingFunctionWithWall --- */

/* Constructor/Destructor */

WeightingFunctionWithWall::WeightingFunctionWithWall(Node *parent): WeightingFunction(parent) {
  init();
}


WeightingFunctionWithWall::~WeightingFunctionWithWall()
{
}


/* Methods */

void WeightingFunctionWithWall::init() {
  m_properties.setClassName("WeightingFunctionWithWall");

  INTPC
    (wallDir, m_wall_dir, -1,
     "Direction in which to find the wall: 0 = x, 1 = y, 2 = z.");

  m_properties.addProperty
    ("leftWall", PropertyList::DOUBLE, &m_left_wall, NULL,
     "Position of the wall to the left.");

  m_properties.addProperty
    ("rightWall", PropertyList::DOUBLE, &m_right_wall, NULL,
     "Position of the wall to the right.");

  m_wall_dir = 0;

  m_left_wall = -10;
  m_right_wall = 10;
}

