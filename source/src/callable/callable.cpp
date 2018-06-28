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



#include "phase.h"
#include "random.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


REGISTER_SMART_ENUM
(Callable_Factory,
 "Callables are generic objects that are called in every simulation "
 "step. Currently (2018-05-04), a Callable is called after the final "
 "integration step within each time step. Refer to the output of "
 "'sympler --help workflow' for more details.\n"
 "The two major differences to Symbols are:\n"
 "  - The computational routine of each Callable handles loops over\n"
 "    particles or pairs of particles by itself, while Symbols are\n"
 "    called *within* a loop.\n"
 "  - Callables do not check for dependencies of their input\n"
 "    variables on output of other Callables. Hence Callables also do\n"
 "    not arrange themselves in corresponding stages to determine\n"
 "    the order in which they will be called. Instead, the order of\n"
 "    execution is strictly the order given as user input. The user\n"
 "    must make sure that the resulting dependencies make sense."
 );


Callable::Callable(Simulation* sim)
  : NodeManyChildren((Node*) sim)
{
}


Callable::~Callable()
{
}


Node* Callable::instantiateChild(const string& name) 
{
  throw gError("Callable::instantiateChild", "Should not be called!");

  return NULL;
}


