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


#include "val_calculator_dirichlet_BC_vector.h"

const SymbolRegister<ValCalculatorDirichletBCVector> val_calc_dirichlet_BC_vector("DirichletBCVector");

using namespace std;

ValCalculatorDirichletBCVector::ValCalculatorDirichletBCVector(string symbol)
  : ValCalculatorDirichletBCArbitrary(symbol)
{}

ValCalculatorDirichletBCVector::ValCalculatorDirichletBCVector(/*Node*/Simulation* parent)
  : ValCalculatorDirichletBCArbitrary(parent)
{
  init();
}

void ValCalculatorDirichletBCVector::init()
{
  m_properties.setClassName("DirichletBCVector");

  m_properties.setDescription(
  		"Saves the pair-specific value of an arbitrary vector of the boundary "
  		"particle used for applying a Dirichlet boundary condition (BC) in each "
  		"pair of particles. This calculator uses a linear approximation as "
  		"described in [70]. The actual value of the Dirichlet boundary "
  		"condition is assumed to be stored in the symbol of the respective "
  		"boundary particle given by the attribute 'variable'.\n"
  		"[70]: J. P. Morris, P. J. Fox, Y. Zhu, J. Comp. Phys. 136 (1997) "
  		"214–226.");

#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorDirichletBCVector::setup()
{
  m_datatype = DataFormat::POINT;
  
  ValCalculatorDirichletBCArbitrary::setup();
}

