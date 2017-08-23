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
  #include <freesteam/region3.h>
  #include <freesteam/b23.h>
}

#include "pressure_calculationTest.h"
#include "pressure_calculation.h"
#include "simulation.h"
#include <vector>

CPPUNIT_TEST_SUITE_REGISTRATION (PressureCalculationTest);

void PressureCalculationTest  :: setUp (void)
{
  /*!
   * Initialize objects
   */
 simulation = new Simulation();
 ps = new PressureCalculation(simulation);
}

void PressureCalculationTest  :: tearDown (void) 
{
  /*!
   * Delete objects
   */
  delete simulation; 
  delete ps;
}
void PressureCalculationTest  :: setupLUTTest (void) {

  double m_Tmin = 650;
  double m_rhomin = 300;
  double m_Tmax= 700;
  double m_rhomax = 600;
  int m_arraysize_density = 100;
  int m_arraysize_temperature = 100;
  ps -> setupLUT(m_Tmin, m_rhomin, m_Tmax, m_rhomax,m_arraysize_density, m_arraysize_temperature);
  double b23Pressure = freesteam_b23_p_T(m_Tmin);
  SteamState S = freesteam_set_pT(b23Pressure, m_Tmin);
  double densityBoundary = freesteam_rho(S);
  if (m_rhomin > densityBoundary) {
    double pressure = freesteam_region3_p_rhoT(m_rhomin , m_Tmin );
    CPPUNIT_ASSERT_EQUAL (ps -> m_array_p[0][0], pressure);
  }
  double m_calcstepRho= (m_rhomax-m_rhomin)/99;
  double m_calcstepT= (m_Tmax-m_Tmin)/99;
  double pressure = freesteam_region3_p_rhoT(600.000000 , 700.000000 );
  CPPUNIT_ASSERT_EQUAL (ps -> m_array_p[99][99], pressure);
}

void PressureCalculationTest  :: calculatePressureTest (void)
{
  double m_Tmin = 650;
  double m_rhomin = 300;
  double m_Tmax= 750;
  double m_rhomax = 600;
  int m_arraysize_density = 1000;
  int m_arraysize_temperature = 1000;
  ps -> setupLUT(m_Tmin, m_rhomin, m_Tmax, m_rhomax,m_arraysize_density, m_arraysize_temperature);
  double pressure = freesteam_region3_p_rhoT(350 , 700 );
  double interpol = ps -> calculatePressure(700, 350, m_Tmin, m_rhomin);
  CPPUNIT_ASSERT_DOUBLES_EQUAL (interpol, pressure,100);
  pressure = freesteam_region3_p_rhoT(600 , 750 );
  interpol = ps -> calculatePressure(750, 600, m_Tmin, m_rhomin);
  CPPUNIT_ASSERT_DOUBLES_EQUAL (interpol, pressure,100);
  pressure = freesteam_region3_p_rhoT(300 , 650 );
  interpol = ps -> calculatePressure(650, 300, m_Tmin, m_rhomin);
  CPPUNIT_ASSERT_DOUBLES_EQUAL (interpol, pressure,100);
}

