/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
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



#ifndef __GRID_METER_SCALAR_TIMES_SCALAR_SUM_H
#define __GRID_METER_SCALAR_TIMES_SCALAR_SUM_H 

#include "grid_meter.h"



/*!
 * Determines the distribution of the user-defined exponent of a scalar quantity
 */
class GridMeterScalarTimesScalarSum: public GridMeter
{
protected:
  /*!
   * Tag offset of the first scalar quantity
   */
  int m_scalar1_o;

  /*!
   * Tag offset of the second scalar quantity
   */
  int m_scalar2_o;

  /*!
   * Offset in the output data for the corresponding result
   */
  int m_result_o;

  //  /*!
  //   * Which exponent to calculate?
  //   */
  //  int m_exponent;

  /*!
   * The name of the first scalar quantity
   */
  string m_scalar1_name;

  /*!
   * The name of the second scalar quantity
   */
  string m_scalar2_name;

 private:

  /*!
   * helper for computing
   */
  vector<double> m_value1;

  /*!
   * helper for computing
   */
  vector<double> m_value2;

  /*!
   * Initialize the property list
   */ 
  void init();

public:
  /*!
   * Constructor
   * @param averager Pointer to the parent \a GridAverager
   */
  GridMeterScalarTimesScalarSum(GridAverager *averager);

  /*!
   * Register the measured quantities in the output data
   */
  virtual void setup();

  /*!
   * Carry out the measurement
   */
  virtual void measure(data_sp data);
};


#endif
