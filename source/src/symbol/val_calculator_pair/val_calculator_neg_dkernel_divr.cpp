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


#include "val_calculator_neg_dkernel_divr.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

const SymbolRegister<ValCalculatorNegDKernelDivr> val_calc_neg_dkernel_divr("ValCalculatorNegDKernelDivr");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorNegDKernelDivr::ValCalculatorNegDKernelDivr(WeightingFunction *wf, string symbol)
  : ValCalculatorPair(symbol), m_wf(wf)
{
//   MSG_DEBUG("ValCalculatorNegDKernelDivr::ValCalculatorNegDKernelDivr", "CONSTRUCTOR");
}

ValCalculatorNegDKernelDivr::ValCalculatorNegDKernelDivr(/*Node*/Simulation* parent)
  : ValCalculatorPair(parent)
{
  m_stage = 0;
  init();
}

void ValCalculatorNegDKernelDivr::init()
{
  m_properties.setClassName("ValCalculatorNegDKernelDivr");

  m_properties.setDescription("Saves -(1/r)*(dW/dr) for each pair of particles, where W is an interpolation function (the kernel) and r is the distance.");
  
  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the computed value.");
  
  STRINGPC
      (weightingFunction, m_wfName,
       "The weighting function to be used.");
  
  m_wfName = "default";
  m_symbolName = "F";

#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorNegDKernelDivr::setup()
{
  m_wf = M_SIMULATION->findWeightingFunction(m_wfName);
  
//  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  ValCalculatorPair::setup();
  
}

void /*pair<size_t, size_t>*/ ValCalculatorNegDKernelDivr::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
{
  m_slot = slot = cp->tagFormat().addAttribute
      ("ValCalculator_" + myName() + "_" + cp->toString(), DataFormat::DOUBLE, false, "F").offset;
}

