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



#ifndef __GRID_AVERAGER_CIRCULAR_H
#define __GRID_AVERAGER_CIRCULAR_H

#include "grid_averager.h"

using namespace std;

/*--- GridAveragerCircular ---*/

/*!
 * For calculating averages in connection with \a BoundaryCouette
 */
class GridAveragerCircular : public GridAverager
{
protected:
  int m_n_z, m_n_r, m_n_phi;
  double m_inner_radius, m_outer_radius;

  double m_diff_r; /* = (m_outer_radius-m_inner_radius)/m_n_cells */

  data_sp m_data;
  vector<double> m_volume; /* Cells have different volumes */

  void findLocations();

  void init();

public:
  GridAveragerCircular(Simulation *simulation);
  virtual ~GridAveragerCircular();

  virtual void measureNow(const double& time);

  virtual void aboutToStart();
	
	virtual void setup();

  virtual void flush();

  virtual double volume(int i) const {
    return m_volume[i];
  }
};

#endif
