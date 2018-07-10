/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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



#include "manager_cell.h"
#include "val_calculator.h"
#include "val_calculator_part.h"
#include "simulation.h"
#include "colour_pair.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorPart::ValCalculatorPart(/*Node*/Simulation* parent): ValCalculator(parent)
{
  init();
}

void ValCalculatorPart::init()
{
  m_properties.setClassName("ValCalculatorPart");
 
  // START: unfinished stuff from 2018-05-08 ////////////////////
  
  // BOOLPC
  //   (selfReset, m_selfReset, "Should this module protect its computed "
  //    "symbol from automatic resetting (to zero) and reset it by "
  //    "itself? The self-reset will be done"
  //    "immediately before the start of computations by any Symbols, "
  //    "including this one in stage 0 and/or 1 as selected. This may be "
  //    "useful to prevent a too early reset, for example by a triggered "
  //    "position update or neighbour list rebuild. In selfReset "
  //    "mode, overwriting existing symbols "
  //    "('overrite = \"yes\"') by this module is forbidden.");

  // m_selfReset = false;

  // END: unfinished stuff from 2018-05-08 ////////////////////

#ifdef _OPENMP
  m_copy_slots.resize(global::n_threads);
#endif

}

void ValCalculatorPart::setup()
{
  ValCalculator::setup();
  
  // START: unfinished stuff from 2018-05-08 ////////////////////
  
  // if(m_selfReset) {
  //   if(m_overwrite)
  //     // Not all children allow to control m_overwrite via XML-input. 
  //     // But since the default in class Symbol is m_overwrite = false,
  //     // this exception and its message should always make sense when
  //     // thrown, since m_overwrite = true is either set via XML or via
  //     // an error in the implementation of a child class. The latter
  //     // should only be metioned in this comment and *not* in the
  //     // exception message, such that we don't confuse users.
  //     throw gError("ValCalculatorPart::setup for module " + className(),
  // 		   "'overwrite = \"yes\"' is not allowed together with "
  // 		   "'selfReset = \"yes\"'.");
  //   m_persistency = true;

  //   // The resetting will be done by an overriden Node::precompute()
  //   if(m_phase == 0) M_CONTROLLER->registerForPrecomputation_0(this);
  //   if(m_phase == 1) M_CONTROLLER->registerForPrecomputation(this);
  //   if(m_phase == 2) {
  //     M_CONTROLLER->registerForPrecomputation(this);
  //     M_CONTROLLER->registerForPrecomputation_0(this);
  //   }
  // } // end of if if(m_selfReset)
  // else
  //   m_persistency = false;

  // END: unfinished stuff from 2018-05-08 ////////////////////

}


