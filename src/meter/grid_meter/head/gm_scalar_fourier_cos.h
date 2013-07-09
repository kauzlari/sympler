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



#ifndef __GRID_METER_SCALAR_FOURIER_COS_H
#define __GRID_METER_SCALAR_FOURIER_COS_H 

#include "gm_scalar.h"



/*!
 * Determines a moment of the distribution of a scalar quantity
 */
class GridMeterScalarFourierCos: public GridMeterScalar
{
protected:

  /*!
   * The wave-vector WITHOUT the factor 2Pi which will automatically be appended
   */
  point_t m_kVec;


#if 0

  /*!
   * The wave length of the plane wave in cartesian direction \a m_dir
   */ 
  double m_lambda;

  /*!
   * Helper: The wave number, i.e., the component of the wave vector in cartesian direction \a m_dir. This is computed from m_lambda
   */ 
  double m_k;

  /*!
   * The cartesian direction of the wave described by the wave number \a m_k
   */ 
  int m_dir;

#endif

  /*!
   * Initialize the property list
   */ 
  void init();

public:
  /*!
   * Constructor
   * @param averager Pointer to the parent \a GridAverager
   */
  GridMeterScalarFourierCos(GridAverager *averager);

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
