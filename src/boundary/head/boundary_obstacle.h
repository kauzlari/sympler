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



#ifndef __BOUNDARY_OBSTACLE_H
#define __BOUNDARY_OBSTACLE_H

#include "cell.h"
#include "boundary_cuboid.h"


/*!
 * A simple cuboid boundary with a cuboid obstacle in the middle
 */
class BoundaryObstacle: public BoundaryCuboid
{
protected:
  /*!
   * Center of the obstacle
   */
  point_t m_obstacle_center;

  /*!
   * Size of the obstacle
   */
  point_t m_obstacle_size;

  /*!
   * Initialize the property list
   */
  void init();
	
public:
  /*!
   * Constructor
   * @param phase Pointer to the \a Phase object this \a Boundary belongs to
   */ 
  BoundaryObstacle(Phase *phase);

  /*!
   * Destructor
   */
  virtual ~BoundaryObstacle();

  /*!
   * Create the boundary in terms of a triangulated surface
   */
  virtual void setup();
};

#endif
