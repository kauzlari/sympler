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
// #include "pressure_calculation.h"
// #include "simulation.h"

CPPUNIT_TEST_SUITE_REGISTRATION (PressureCalculationTest);

void PressureCalculationTest  :: setUp (void)
{
  /*!
   * Initialize objects
   */
  m_simulation = new Simulation();
  m_ps = new FakePressureCalculation(m_simulation);

  // Minimum and maximum Inputvalues
  m_Tmin = 650;
  m_rhomin = 300;
  m_Tmax= 750;
  m_rhomax = 600;
  // Size of the arrays
  m_arraysize_density = 10;
  m_arraysize_temperature = 10;
  
  m_ps->set_m_Tmax(m_Tmax);
  m_ps->set_m_rhomax(m_rhomax);
  m_ps->set_m_Tmin(m_Tmin);
  m_ps->set_m_rhomin(m_rhomin);
  m_ps->set_m_arraysize_temperature(m_arraysize_temperature);
  m_ps->set_m_arraysize_density(m_arraysize_density);
  
  m_ps ->setupLUT();

}

void PressureCalculationTest  :: tearDown (void) 
{
  /*!
   * Delete objects
   */
  delete m_simulation; 
  delete m_ps;
}

void PressureCalculationTest  :: setupLUTTest (void) {
  
  // Initialization of the LUT now done in setUp of test class
  // m_ps -> setupLUT();

  double b23Pressure = freesteam_b23_p_T(m_Tmin);
  SteamState S = freesteam_set_pT(b23Pressure, m_Tmin);
  double densityBoundary = freesteam_rho(S);
  if (m_rhomin > densityBoundary) {
    double pressure = freesteam_region3_p_rhoT(m_rhomin, m_Tmin);
    // Assertions to check if the first and the last entry of the array are correct.
    CPPUNIT_ASSERT_EQUAL (m_ps -> returnLUTvals()[0][0], pressure);
  }
  else
    // case should not happen, so check for it
    CPPUNIT_ASSERT_EQUAL (2.0, 1.0);
    
  double pressure = freesteam_region3_p_rhoT(m_rhomax, m_Tmax);

  CPPUNIT_ASSERT_EQUAL (m_ps -> returnLUTvals()[m_arraysize_density - 1][m_arraysize_temperature - 1], pressure);

}

void PressureCalculationTest  :: calculatePressureTest (void)
{
  execPressureTest(350., 700.);
  execPressureTest(600., 750.);
  execPressureTest(300., 650.); 
}

void PressureCalculationTest  :: execPressureTest(const double& density, const double& temperature) {
  
  // Assertion to check if the interpolation is correct.
  double interpol;
  double pressure;
  string exceptMsg;
  
  exceptMsg = "ERROR: Unspecified exception in PressureCalculationTest::calculatePressureTest. Please check!";
  pressure = freesteam_region3_p_rhoT(350., 700.);
  try {
    interpol = m_ps -> calculatePressure(700., 350.);
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }
  CPPUNIT_ASSERT_NO_THROW_MESSAGE(exceptMsg, m_ps -> calculatePressure(700., 350.));
  CPPUNIT_ASSERT_DOUBLES_EQUAL (interpol, pressure, interpol/100.); // relativ delta: 1% of expected
  
}