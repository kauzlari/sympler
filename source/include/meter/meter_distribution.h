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



// should be usable for all dimensionalities
#ifndef __METER_DISTRIBUTION_H
#define __METER_DISTRIBUTION_H

#include "meter.h"

using namespace std;

//---- newclass -----------------------------------------------------------

class Simulation;


/*!
 * Determine the one particle distribution function of a quantity
 */
class MeterDistribution : public Meter
{
protected:
  /*!
   * The species for which to perform measurements
   */
  string m_species;

  /*!
   * The color for which to perform measurements
   */
  size_t m_colour;

  /*!
   * The number of bins to use for the distribution
   */
  int m_nbins;

  /*!
   * The number of particle that were binned
   */
  size_t m_n_particles;

  /*!
   * Tag offset of the data for which to calculate the distribution
   */
  int m_offset;

  /*!
   * Current simulation step
   */
  int m_step;

  /*!
   * Average over this number of steps
   */
  int m_avg_over;
  
  /*!
   * Value corresponding to the left bin
   */
  double m_min;

  /*!
   * Value corresponding to the right bin
   */
  double m_max;

  /*!
   * = \a m_max - \a m_min
   */
  double m_size;

  /*!
   * Spacing of the bins
   */
  double m_bin_spacing;

  /*!
   * What to measure
   */
  string m_what;

  /*!
   * Values of the bins
   */
  vector_double_sp m_bins;

  /*!
   * Positions of the bins
   */
  vector_double_sp m_bin_positions;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  MeterDistribution(Simulation *simulation);

  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   * @param everyN Measure only every other time step
   */
  MeterDistribution(Simulation* simulation, const size_t& everyN/*, bool only*/);

  /*!
   * Destructor
   */
  virtual ~MeterDistribution();

  /*!
   * Calculate the total number of cells
   */
  virtual void setup();

  /*!
   * Determine the averages for this timestep
   */
  virtual void measureNow(const double& time);
};

#endif
