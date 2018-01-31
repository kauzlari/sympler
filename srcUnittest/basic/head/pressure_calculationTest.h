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



#ifndef PRESSURECALCULATIONTEST_H
#define PRESSURECALCULATIONTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "pressure_calculation.h"
#include "simulation.h"
using namespace std;



/*!
 * Fake class for testing class \a PressureCalculation. Current resons 
 * for its need:
 * - private/protected member m_Tmax, m_pmax should be possible to set
 *   by public method
 */
class FakePressureCalculation : public PressureCalculation
{
  
 public: 
  
  FakePressureCalculation (Simulation* parent)
    : PressureCalculation (parent)
  {}
  
  /*!
   * Destructor
   */
  virtual ~FakePressureCalculation() {
  }
  
  /*!
   * Set private member \a m_Tmax to given value
   * @param Tmax value that m_Tmax should be set to
   */
  void set_m_Tmax(const double& Tmax) {
    m_Tmax = Tmax;
  }
  
  /*!
   * Set private member \a m_rhomax to given value
   * @param rhomax value that m_rhomax should be set to
   */
  void set_m_rhomax(const double& rhomax) {
    m_rhomax = rhomax;
  }
  
  /*!
   * Set private member \a m_Tmin to given value
   * @param Tmin value that m_Tmin should be set to
   */
  void set_m_Tmin(const double& Tmin) {
    m_Tmin = Tmin;
  }
  
  /*!
   * Set private member \a m_rhomin to given value
   * @param rhomin value that m_rhomin should be set to
   */
  void set_m_rhomin(const double& rhomin) {
    m_rhomin = rhomin;
  }
  
  /*!
   * Set private member \a m_arraysize_temperature to given value
   * @param size value that m_arraysize_temperature should be set to
   */
  void set_m_arraysize_temperature(size_t size) {
    m_arraysize_temperature = size;
  }
  
  /*!
   * Set private member \a m_arraysize_pressure to given value
   * @param size value that m_arraysize_pressure should be set to
   */
  void set_m_arraysize_density(size_t size) {
    m_arraysize_density = size;
  }
    
}; // end of FakePressureCalculation


class PressureCalculationTest : public CPPUNIT_NS :: TestFixture
{
  CPPUNIT_TEST_SUITE (PressureCalculationTest);
  CPPUNIT_TEST (setupLUTTest);
  CPPUNIT_TEST (calculatePressureTest);
  CPPUNIT_TEST_SUITE_END ();
  
 public:

  void setUp (void);
  void tearDown (void);

 protected:

  void setupLUTTest (void);
  void calculatePressureTest (void);

 private:

  Simulation *m_simulation; 
  FakePressureCalculation *m_ps;

  double m_Tmax;
  double m_rhomax;
  double m_Tmin;
  double m_rhomin;
  size_t m_arraysize_temperature;
  size_t m_arraysize_density;

  /*!
   * Helper for multiple execution
   * @param density Density of the thermodynamic state
   * @param tempersture Temperature of the thermodynamic state
   */
  void execPressureTest(const double& density, const double& temperature);
  
};

#endif
