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



#ifndef __GRID_METER_H
#define __GRID_METER_H

#include "node.h"
#include "smart_enum.h"
#include "data_format.h"


/*--- GridGridMeter ---*/

class GridAverager;

#define M_GRID_AVERAGER    ((GridAverager*) m_parent)


/*!
 * Base class for all \a GridMeter s that measure quantities 
 * at a certain position in the domain. They are used together with a \a GridAverager
 * which defines the grid for which to determine quantities
 */
class GridMeter: public Node
{
protected:
  /*!
   * Species to measure
   */
  string m_species;

  /*!
   * Color to measure
   */
  size_t m_colour;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param averager Pointer to the parent \a GridAverager
   */
  GridMeter(GridAverager *averager);

  /*!
   * Destructor
   */
  virtual ~GridMeter();

  /*!
   * Get the color corresponding to the given species
   */
  virtual void setup();

  /*!
   * Make the measurement
   */
  virtual void measure(data_sp data) = 0;

  /*!
   * This averaging step has finished
    * @param n_steps is needed for correct computation of variances
   */
  virtual void finishStep(data_sp data, size_t n_steps) const {
  }


  /* --- Utility functions ---*/

  /*!
   * Calculate the variance from the mean and the mean squared values
   * @param values The mean values
   * @param values_sq The mean squared values
   * @param n_steps is needed for correct computation of variances
   */
  static void calcVariance(const vector_double_sp &values, vector_double_sp &values_sq, size_t n_steps);

  /*!
   * Calculate the variance from the mean and the mean squared values
   * @param values The mean values
   * @param values_sq The mean squared values
   * @param n_steps is needed for correct computation of variances
   */
  static void calcVariance(const vector_point_sp &values, vector_point_sp &values_sq, size_t n_steps);
};



/*--- Factories ---*/

class GridMeter_Factory: public SmartEnum<GridMeter_Factory>
{
public:
    virtual GridMeter *instantiate(GridAverager *averager) const = 0;

protected:
    GridMeter_Factory(const string &name)
    : SmartEnum<GridMeter_Factory>(name) { }
};


template <class T>
class GridMeter_Register: public GridMeter_Factory
{
public:
    GridMeter_Register(const string &name)
    : GridMeter_Factory(name) { }

    virtual GridMeter *instantiate(GridAverager *averager) const;
};



//---- Inline functions ----

template <class T>
inline GridMeter *GridMeter_Register<T>::instantiate(GridAverager *averager) const
{
    return new T(averager);
}

#endif
