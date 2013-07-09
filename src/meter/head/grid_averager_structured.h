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



#ifndef __GRID_AVERAGER_STRUCTURED_H
#define __GRID_AVERAGER_STRUCTURED_H

#include "grid_averager.h"

using namespace std;

/*!
 * Create an evenly spaced (structured) grid
 */
class GridAveragerStructured : public GridAverager
{
 protected:
  /*!
   * Number of cells in x, y and z direction
   */
  int_point_t m_n;

  /*!
   * The actual data
   */
  data_sp m_data;

  /*!
   * The volume of a cell. All cells have the same volume.
   */
  double m_volume;

  /*!
   * Initialize the location array
   */
  void findLocations();

  /*!
   * Initialize the property list
   */
  void init();

 public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  GridAveragerStructured(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~GridAveragerStructured();

  /*!
   * Determine the averages for this timestep
   */
  virtual void measureNow(const double& time);

  /*!
   * Initialize the cell volume \a m_volume
   */
  virtual void aboutToStart();
	
  /*!
   * Calculate the total number of cells
   */
  virtual void setup();

  /*!
   * Flush the information at the end of the simulation
   */
  virtual void flush();

  /*!
   * Return the volume of cell \a i
   * @param i Index of the cell
   */
  virtual double volume(int i) const {
    return m_volume;
  }
};

#endif
