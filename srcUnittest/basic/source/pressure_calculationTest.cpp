/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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

extern "C" {
  #include <freesteam/steam_pT.h>
  #include <freesteam/region3.h>
  #include <freesteam/b23.h>
}

#include "pressure_calculationTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION (PressureCalculationTest);

void PressureCalculationTest  :: setUp (void)
{
  /*!
   * Initialize objects
   */
  m_simulation = new Simulation();
  m_pc = new PressureCalculation(m_simulation);

  // Minimum and maximum Inputvalues
  // density
  m_var1Min = 300;
  m_var1Max = 600;
  // temperature
  m_var2Min = 650;
  m_var2Max= 750;

  PCacheIAPWSIF97Test::setUp();  
}

void PressureCalculationTest  :: tearDown (void) 
{
  PCacheIAPWSIF97Test::tearDown();
}

void PressureCalculationTest  :: calculateResultTest (void)
{
  // 1st arg: density, 2nd arg: temperature 
  execCalculationTest(350., 700.);
  execCalculationTest(600., 750.);
  execCalculationTest(300., 650.); 
}
