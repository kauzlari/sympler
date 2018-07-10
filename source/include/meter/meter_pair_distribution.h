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



#ifndef __METER_PAIR_DISTRIBUTION_H
#define __METER_PAIR_DISTRIBUTION_H

#include "meter.h"
#include "colour_pair.h"

using namespace std;

//---- newclass -----------------------------------------------------------

class Simulation;

/*!
 * Determine the pair distribution function g(r)
 */
class MeterPairDistribution : public Meter
{
protected:
  /*!
   * Species for which to determine the pair distribution
   */
  string m_species;

  /*!
   * The corresponding color
   */
  size_t m_colour;

  /*!
   * The corresponding color pair (where both colors are actually the same)
   */
  ColourPair *m_cp;

  /*!
   * Number of bins to use
   */
  int m_nbins;
  
  /*!
   * The current simulation step
   */
  int m_step;

  /*!
   * Average over this many of simulation step
   */
  int m_avg_over;

  /*!
   * The maximum radius taken for the pair distribution function
   */
  double m_max_radius;

  /*!
   * The bin size
   */
  double m_bin_size;

  /*!
   * The reciprocal of the bin size
   */
  double m_r_bin_size;

  /*!
   * The measured values
   */
  vector_double_sp m_bins;

  /*!
   * The variances of the values
   */
  vector_double_sp m_variances;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  MeterPairDistribution(Simulation *simulation);

  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   * @param everyN Measure only every other time step
   */
  MeterPairDistribution(Simulation* simulation, const size_t& everyN/*, bool only*/);

  /*!
   * Destructor
   */
  virtual ~MeterPairDistribution();

  /*!
   * Get the \a ColourPair, etc.
   */
  virtual void setup();

  /*!
   * Make the measurement
   */
  virtual void measureNow(const double& time);
};

#endif
