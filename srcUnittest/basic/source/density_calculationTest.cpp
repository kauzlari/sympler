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
#include "density_calculationTest.h"
// #include "density_calculation.h"
// #include "simulation.h"

CPPUNIT_TEST_SUITE_REGISTRATION (DensityCalculationTest);

void DensityCalculationTest  :: setUp (void) {
  /*!
   * Initialize objects
   */
 m_simulation = new Simulation();
 m_ps = new FakeDensityCalculation(m_simulation);

 m_Tmax = 700.;
 m_pmax = 25000000.;
 m_Tmin = 550.;
 m_pmin = 23000000.;
 m_arraysize_temperature = 10;
 m_arraysize_pressure = 10;

 m_ps->set_m_Tmax(m_Tmax);
 m_ps->set_m_pmax(m_pmax);
 m_ps->set_m_Tmin(m_Tmin);
 m_ps->set_m_pmin(m_pmin);
 m_ps->set_m_arraysize_temperature(m_arraysize_temperature);
 m_ps->set_m_arraysize_pressure(m_arraysize_pressure);

 // FIXME: does not need any arguments!!!
 m_ps ->setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax, m_arraysize_pressure, m_arraysize_temperature);
 
}

void DensityCalculationTest  :: tearDown (void) {
  /*!
   * Delete objects
   */
  delete m_simulation;
  delete m_ps;
}

void DensityCalculationTest  :: setupLUTTest (void) {

  // Initialization of the LUT
  // m_ps ->setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax,m_arraysize_pressure, m_arraysize_temperature);
  SteamState S = freesteam_set_pT(m_pmin, m_Tmin);
  double density = freesteam_rho(S);
  // Assertions to check if the first and the last content of the array are correct.
  CPPUNIT_ASSERT_EQUAL (m_ps -> m_array_rho[0][0], density);
  // Calculation steps between expansion points
  double m_calcstepP= (m_pmax-m_pmin)/(m_arraysize_pressure-1);
  double m_calcstepT= (m_Tmax-m_Tmin)/(m_arraysize_temperature-1);
  S = freesteam_set_pT(m_pmin+ (m_arraysize_pressure-1)*m_calcstepP, m_Tmin+ (m_arraysize_temperature-1)*m_calcstepT);
  density = freesteam_rho(S);
  CPPUNIT_ASSERT_EQUAL (m_ps -> m_array_rho[(m_arraysize_pressure-1)][(m_arraysize_temperature-1)], density);
}

void DensityCalculationTest  :: calculateDensityTest (void)
{

  // Initialization of the LUT
  m_ps -> setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax, m_arraysize_pressure, m_arraysize_temperature );
  SteamState S = freesteam_set_pT(24000000, 650);
  double density = freesteam_rho(S);
  // Assertion to check if the interpolation is correct.
  double temp;
  string tempString = "ERROR: Unspecified bug in DensityCalculationTest::calculateDensityTest. Please check!";
  gError tempGerror;

  try {
    temp = m_ps -> calculateDensity(650, 24000000, m_Tmin,m_pmin);
  } catch(gError& err) {

  	tempString = err.message(); 

  }

  CPPUNIT_ASSERT_NO_THROW_MESSAGE(tempString, m_ps -> calculateDensity(650, 24000000, m_Tmin,m_pmin));

  CPPUNIT_ASSERT_DOUBLES_EQUAL (temp, density, 0.2);
  // CPPUNIT_ASSERT_DOUBLES_EQUAL (ps -> calculateDensity(650, 24000000, m_Tmin,m_pmin), density, 0.1);

}


