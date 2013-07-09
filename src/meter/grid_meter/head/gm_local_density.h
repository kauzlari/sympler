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



#ifndef __GRID_METER_LOCAL_DENSITY_H
#define __GRID_METER_LOCAL_DENSITY_H 

#include "grid_meter.h"
// #include "colour_pair.h"
#include "manager_cell.h"
#include "weighting_function.h"


/*--- GridMeterLocalDensity ---*/

/*!
 * Determines the average local density per particle
 */
class GridMeterLocalDensity: public GridMeter
{
protected:
//  /*!
//   * The color pair for which to determine the local density
  //   */
//  ColourPair *m_cp;

  /*!
   * Tag offset for the local density
   */
  size_t m_density_o;

  /*!
   * Tag offset in the output data for the averaged local density
   */
  size_t m_density_m_o;

  /*!
   * Tag offset in the output data for the variance of the local density
   */
  size_t m_density_v_o;

//  /*!
//   * The name of the weighting function to use for the local density calculation
  //   */
//  string m_weighting_function;

//  /*!
//   * A pointer to the actual weighting function impelementation
  //   */
//  WeightingFunction *m_wf;

  /*!
  * Density to be measured
  */
  string m_rhoSymbol;
  
  /*!
   * Will the local density be computed over all ColourPairs (CP) or only over the CP 
   * corresponding to the chosen species?
   */
  bool m_oneProp;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param averager Pointer to the parent \a GridAverager
   */
  GridMeterLocalDensity(GridAverager *averager);

  /*!
   * Register the measured quantities in the output data
   */
  virtual void setup();

  /*!
   * Carry out the measurement
   */
  virtual void measure(data_sp data);

  /*!
   * Calculate the variances
   */
  virtual void finishStep(data_sp data, size_t n_steps) const {
    calcVariance(data->vectorDoubleByOffset(m_density_m_o), data->vectorDoubleByOffset(m_density_v_o), n_steps);
  }
};


#endif
