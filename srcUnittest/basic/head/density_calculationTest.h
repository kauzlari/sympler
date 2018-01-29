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



#ifndef DENSITYCALCULATIONTEST_H
#define DENSITYCALCULATIONTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "density_calculation.h"
#include "simulation.h"
using namespace std;


/*!
 * Fake class for testing class \a Densityclaculation. Current resons 
 * for its need:
 * - private/protected member m_Tmax, m_pmax should be possible to set
 *   by public method
 */
class FakeDensityCalculation : public DensityCalculation
{
  
 public: 
  
  FakeDensityCalculation (Simulation* parent)
    : DensityCalculation (parent)
  {}
  
  /*!
   * Destructor
   */
  virtual ~FakeDensityCalculation() {
  }
  
  /*!
   * Set private member \a m_Tmax to given value
   * @param Tmax value that m_Tmax should be set to
   */
  void set_m_Tmax(const double& Tmax) {
    m_Tmax = Tmax;
  }
  
  /*!
   * Set private member \a m_pmax to given value
   * @param pmax value that m_pmax should be set to
   */
  void set_m_pmax(const double& pmax) {
    m_pmax = pmax;
  }
  
  /*!
   * Set private member \a m_Tmin to given value
   * @param Tmin value that m_Tmin should be set to
   */
  void set_m_Tmin(const double& Tmin) {
    m_Tmin = Tmin;
  }
  
  /*!
   * Set private member \a m_pmin to given value
   * @param pmin value that m_pmin should be set to
   */
  void set_m_pmin(const double& pmin) {
    m_pmin = pmin;
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
  void set_m_arraysize_pressure(size_t size) {
    m_arraysize_pressure = size;
  }
    
};

class DensityCalculationTest : public CPPUNIT_NS :: TestFixture
{
  CPPUNIT_TEST_SUITE (DensityCalculationTest);
  CPPUNIT_TEST (setupLUTTest);
  CPPUNIT_TEST (calculateDensityTest);
  CPPUNIT_TEST_SUITE_END ();
  
 public:
  
  void setUp (void);
  void tearDown (void);
  
 protected:
  
  void setupLUTTest (void);
  void calculateDensityTest (void);
  
 private:
  
  Simulation *m_simulation; 
  FakeDensityCalculation *m_ps;

  double m_Tmax;
  double m_pmax;
  double m_Tmin;
  double m_pmin;
  size_t m_arraysize_temperature;
  size_t m_arraysize_pressure;

};
#endif
