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



#ifndef __GRID_METER_SELF_DIFFUSION_H
#define __GRID_METER_SELF_DIFFUISON_H

#include "grid_meter.h"



/*--- GridMeterSelfDiffusion ---*/

/*!
 * Determine the self diffusion coefficient.
 * See: D. C. Rapaport, The Art of Molecular Dynamics Simulation, p. 116
 */
class GridMeterSelfDiffusion: public GridMeter
{
protected:
  /*!
   * Tag offset of the displacement
   */
  int m_dispOffset;

  /*!
   * The symbol name of the displacement
   */
  string m_dispName;


  /*!
   * The offset for the self diffusion coefficient
   */
  int m_D_o;

  /*!
   * Is this the first measurement? Needed to know in order to
   * set \a m_initDisp
   */
  bool m_first_measurement;

  /*!
   * Internal helper: initial displacement to be subtracted from measurement; 
   */
  point_t m_initDisp;

  /*!
   * Initialize (i.e., set the property list information)
   */
  void init();

public:
  /*!
   * Constructor
   */
  GridMeterSelfDiffusion(GridAverager *averager);

  /*!
   * Setup before simulation starts
   */
  virtual void setup();
  
  /*!
   * Take a measurement
   */
  virtual void measure(data_sp data);

  /*!
   * Finish this measurement (averaging) step
   */
  virtual void finishStep(data_sp data, size_t n_steps) const {
  }
};


#endif
