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




#ifndef __GRID_METER_VELOCITY_MOMENT_H
#define __GRID_METER_VELOCITY_MOMENT_H

#include "grid_meter.h"



/*--- GridMeterDensity ---*/


/*!
 * Determines a certain moment of the velocity distribution within a cell
 */
class GridMeterVelocityMoment: public GridMeter
{
 protected:
  /*!
   * Order of the moment to determine, i.e. <v> -> order = 1, <v^2> -> order = 2 ...
   */
  int m_order;

  /*!
   * Offset of the moment of the absolute value of the velocity in the output data
   */
  size_t m_offset;

  /*!
   * Offset of the moment of the different velocity components in the output data
   */
  size_t m_offset_by_dir;

  /*!
   * Initialize the property list
   */
  void init();

 public:
  /*!
   * Constructor
   * @param averager Pointer to the parent \a GridAverager
   */
  GridMeterVelocityMoment(GridAverager *averager);

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
