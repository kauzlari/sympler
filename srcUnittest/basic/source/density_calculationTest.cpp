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
}
#include "density_calculationTest.h"
#include "density_calculation.h"
#include "simulation.h"
#include <vector>


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
  double m_Tmin = 500;
  double m_pmin = 23000000;
  double m_Tmax= 700;
  double m_pmax = 24000000;
  int m_arraysize_pressure = 100;
  int m_arraysize_temperature = 100;
  ps ->setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax,m_arraysize_pressure, m_arraysize_temperature);
  SteamState S = freesteam_set_pT(m_pmin, m_Tmin);
  double density = freesteam_rho(S);
  CPPUNIT_ASSERT_EQUAL (ps -> m_array_rho[0][0], density);
  double m_calcstepP= (m_pmax-m_pmin)/99;
  double m_calcstepT= (m_Tmax-m_Tmin)/99;
  S = freesteam_set_pT(m_pmin+ 99*m_calcstepP, m_Tmin+ 99*m_calcstepT);
  density = freesteam_rho(S);
  CPPUNIT_ASSERT_EQUAL (ps -> m_array_rho[99][99], density);
}

void DensityCalculationTest  :: calculateDensityTest (void)
{
  double m_Tmin = 500;
  double m_pmin = 23000000;
  double m_Tmax= 700;
  double m_pmax = 24000000;
  int m_arraysize_pressure = 100;
  int m_arraysize_temperature = 100;
  double m_calcstepP= (m_pmax-m_pmin)/99;
  double m_calcstepT= (m_Tmax-m_Tmin)/99;
  ps -> setupLUT(m_Tmin, m_pmin, m_Tmax, m_pmax, m_arraysize_pressure, m_arraysize_temperature );
  //double ipDensity= ps -> calculateDensity(662, 24000000, m_Tmin, m_pmin);
  SteamState S = freesteam_set_pT(24000000, 700);
  double density = freesteam_rho(S);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ((ps -> calculateDensity(700,24000000 , m_Tmin,m_pmin)), density,0.1);

}


