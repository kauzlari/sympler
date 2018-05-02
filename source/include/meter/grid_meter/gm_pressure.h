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



#ifndef __GRID_METER_PRESSURE_H
#define __GRID_METER_PRESSURE_H

#include "grid_meter.h"



/*--- GridMeterPressure ---*/


/*!
 * Determines the local pressure
 */
class GridMeterPressure: public GridMeter
{
protected:
  /*!
   * Offset of the pressure in the output data
   */
  size_t m_p_offset;
  
  /*!
   * The name of the stress quantity
   */
  string m_stress_name;
  
  /*!
   *  Offset of the stress
   */
  int m_stress_offset;

  /*!
   * Initialize the property list
   */
  void init();

 public:
  /*!
   * Constructor
   * @param averager Pointer to the parent \a GridAverager
   */
  GridMeterPressure(GridAverager *averager);

  /*!
   * Register the measured quantities in the output data
   */
  virtual void setup();

  /*!
   * Carry out the measurement
   */
  virtual void measure(data_sp data);

  /*!
   * Fixme!!! Does nothing right now... Does the variance of a pressure make sense?
   */
  virtual void finishStep(data_sp data, size_t n_steps) const {
    //    calcVariance(data->vectorDoubleByOffset(m_p), data->vectorDoubleByOffset(m_pv));
  }
};


#endif
