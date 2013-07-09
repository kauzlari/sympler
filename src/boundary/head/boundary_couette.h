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



#ifndef __BOUNDARY_COUETTE_H
#define __BOUNDARY_COUETTE_H

#include "cell.h"
#include "boundary_arbitrary.h"

/*!
 * Simple Couette flow geometry
 */
class BoundaryCouette: public BoundaryArbitrary
{
protected:
  /*!
   * Inner tube radius
   */
  double m_inner_radius;

  /*!
   * Outer tube radius
   */
  double m_outer_radius;

  /*!
   * Tube height
   */
  double m_height;

  /*!
   * Is the geometry periodic in direction of the tube?
   */
  bool m_periodic;

  /*!
   * Name of the reflector for the inner tube wall
   */
  string m_inner_reflector;

  /*!
   * Name of the reflector for the outer tube wall
   */
  string m_outer_reflector;

  /*!
   * Number of edge elements to use for the discretization of the curved
   * tube surface.
   */
  int m_n_edge_elements;

  /*!
   * Initialize the property list
   */
  void init();
	
public:
  /*!
   * Constructor
   * @param phase Pointer to the parent \a Phase this \a Boundary belongs to
   */
  BoundaryCouette(Phase *phase);
  
  /*!
   * Destructor
   */
  virtual ~BoundaryCouette();
	
  /*!
   * Initialize cell subdivision
   */
  virtual void setup(Simulation* sim, ManagerCell *mgr);

  /*!
   * Create the geometry in terms of a triangulated surface
   */
  virtual void setup();
};

#endif
