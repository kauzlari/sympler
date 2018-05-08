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

#ifdef HAVE_FREESTEAM

extern "C" {
  #include <freesteam/steam_pT.h>
}

#include "pca_iapws-if97_rhoTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION (PCacheIAPWSIF97rhoTest);

void PCacheIAPWSIF97rhoTest  :: setUp (void)
{
  /*!
   * Initialize objects
   */
  m_simulation = new Simulation();
  m_pc = new PCacheIAPWSIF97rho(m_simulation);

  // Minimum and maximum Inputvalues
  // pressure
  m_var1Min = 23000000.;
  m_var1Max = 25000000.;
  // temperature
  m_var2Min = 550.;
  m_var2Max= 700;

  PCacheIAPWSIF97TwoVarTest::setUp();  
}

void PCacheIAPWSIF97rhoTest  :: tearDown (void) 
{
  PCacheIAPWSIF97TwoVarTest::tearDown();
}

void PCacheIAPWSIF97rhoTest  :: calculateResultTest (void)
{
  // 1st arg: pressure, 2nd arg: temperature 
  execCalculationTest(24000000., 650.);
}

#endif    // HAVE_FREESTEAM
