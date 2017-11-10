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


extern "C" {
  #include <freesteam/steam_pT.h>
}
#include "density_calculationTest.h"
#include "density_calculation.h"
#include "simulation.h"

CPPUNIT_TEST_SUITE_REGISTRATION (DensityCalculationTest);

void DensityCalculationTest  :: setUp (void) {
  /*!
   * Initialize objects
   */
 simulation = new Simulation();
 ps = new DensityCalculation(simulation);
}

void DensityCalculationTest  :: tearDown (void) {
  /*!
   * Delete objects
   */
  delete simulation;
  delete ps;
}

void DensityCalculationTest  :: setupLUTTest (void) {
  // Minimum and maximum Inputvalues
  double m_Tmin = 500;
  double m_pmin = 23000000;
  double m_Tmax= 700;
  double m_pmax = 25000000;
  // Size of the arrays
  int m_arraysize_pressure = 10;
  int m_arraysize_temperature = 10;
  // Initialization of the LUT
  ps ->setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax,m_arraysize_pressure, m_arraysize_temperature);
  SteamState S = freesteam_set_pT(m_pmin, m_Tmin);
  double density = freesteam_rho(S);
  // Assertions to check if the first and the last content of the array are correct.
  CPPUNIT_ASSERT_EQUAL (ps -> m_array_rho[0][0], density);
  // Calculation steps between expansion points
  double m_calcstepP= (m_pmax-m_pmin)/(m_arraysize_pressure-1);
  double m_calcstepT= (m_Tmax-m_Tmin)/(m_arraysize_temperature-1);
  S = freesteam_set_pT(m_pmin+ (m_arraysize_pressure-1)*m_calcstepP, m_Tmin+ (m_arraysize_temperature-1)*m_calcstepT);
  density = freesteam_rho(S);
  CPPUNIT_ASSERT_EQUAL (ps -> m_array_rho[(m_arraysize_pressure-1)][(m_arraysize_temperature-1)], density);
}

void DensityCalculationTest  :: calculateDensityTest (void)
{
  // Minimum and maximum Inputvalues
  double m_Tmin = 550;
  double m_pmin = 23000000;
  double m_Tmax= 700;
  double m_pmax = 24000000;
  // Size of the arrays
  int m_arraysize_pressure = 10;
  int m_arraysize_temperature = 10;
  // Initialization of the LUT
  ps -> setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax, m_arraysize_pressure, m_arraysize_temperature );
  SteamState S = freesteam_set_pT(24000000, 650);
  double density = freesteam_rho(S);
  // Assertion to check if the interpolation is correct.
  CPPUNIT_ASSERT_DOUBLES_EQUAL ((ps -> calculateDensity(650, 24000000, m_Tmin,m_pmin)), density, 0.1);

}


