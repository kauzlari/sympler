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
}

#include "pca_iapws-if97_cpTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION (PCacheIAPWSIF97CpTest);

void PCacheIAPWSIF97CpTest  :: setUp (void)
{
  /*!
   * Initialize objects
   */
  m_simulation = new Simulation();
  m_pc = new PCacheIAPWSIF97Cp(m_simulation);

  // Minimum and maximum Inputvalues
  // temperature
  m_var1Min = 550.;
  m_var1Max= 700.;

  // Arguments must have been set in setUp of subclass before
  // Cast (PCacheIAPWSIF97Cp*) required as long as
  // PCacheIAPWSIF97WithConstTestGetter is only implemented for the
  // non-generalised class PCacheIAPWSIF97Cp
  m_pcWithConstGetter = new PCacheIAPWSIF97WithConstTestGetter((PCacheIAPWSIF97Cp*) m_pc);
  
  // The constant pressure for Cp
  m_pcWithConstGetter->set_m_constP(25000000.);
  
  PCacheIAPWSIF97OneVarTest::setUp();  
}

void PCacheIAPWSIF97CpTest  :: tearDown (void) 
{
  PCacheIAPWSIF97OneVarTest::tearDown();
}

void PCacheIAPWSIF97CpTest  :: calculateResultTest (void)
{

  // 1st argument: temperature 
  execCalculationTest(650.);
}

