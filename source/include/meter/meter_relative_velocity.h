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
#ifndef __METER_RELATIVE_VELOCITY_H
#define __METER_RELATIVE_VELOCITY_H

#include "meter.h"
#include "colour_pair.h"

using namespace std;

//---- Classes ----

class Simulation;

/*!
 * Measures the first four cummulants of the eij*vij distribution
 */
class MeterRelativeVelocity : public Meter
{
protected:
  /*!
   * Species to make the measurements for
   */
  pair<string, string> m_species;

  /*!
   * \a ColourPair corresponding to the species
   */
  ColourPair *m_cp;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  MeterRelativeVelocity(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~MeterRelativeVelocity();

  /*!
   * Gets the \a ColourPair
   */
  virtual void setup();

  /*!
   * Do the measurement
   */
  virtual void measureNow(const double& time);
};

#endif
