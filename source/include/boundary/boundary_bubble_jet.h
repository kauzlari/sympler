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



#ifndef __BOUNDARY_BUBBLE_JET_H
#define __BOUNDARY_BUBBLE_JET_H

#include "cell.h"
#include "boundary_arbitrary.h"

/*!
 * Simple boundary mimicing a bubble jet print head
 */
class BoundaryBubbleJet: public BoundaryArbitrary
{
protected:
  /*!
   * The width of the reservoir
   */
  double m_reservoir_width;

  /*!
   * The width of the inlet
   */
  double m_inlet_width;

  /*!
   * The width of the diffusor
   */
  double m_diffusor_width;

  /*!
   * The height of the reservoir
   */
  double m_reservoir_height;

  /*!
   * The height of the inlet
   */
  double m_inlet_height;

  /*!
   * The height of the diffusor
   */
  double m_diffusor_height;
 
  /*!
   * The axis at which to align the structure
   */
  int m_axis_dir;

  /*!
   * Initialize the property list
   */
  void init();
	
public:
  /*!
   * Constructor
   * @param phase Pointer to the \a Phase object this \a Boundary belongs to
   */
  BoundaryBubbleJet(Phase *phase);
  
  /*!
   * Destructor
   */ 
  virtual ~BoundaryBubbleJet();
	
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

