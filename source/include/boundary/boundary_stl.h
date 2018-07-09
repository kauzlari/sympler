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



#ifndef __BOUNDARY_STL_H
#define __BOUNDARY_STL_H

#include "cell.h"
#include "boundary_with_inlet.h"

/*!
 * Boundary that loads the surface description from a STL file
 */
class BoundarySTL: public BoundaryWithInlet
{
protected:
  /*!
   * STL file name
   */
  string m_filename;

  /*!
   * Invert the normal vector of the STL. Usually needs to be set
   * to true, I never encountered an STL where the normal vectors
   * were pointing to the inside of the geometry, but I wouldn't count 
   * on it given the simplicity/crappyness of the STL format.
   */
  bool m_invert_normals;

  /*!
   * Gives the direction where the inlet can be found: 'x', 'y' or 'z', or 'none'
   */
  string m_inlet_str;

  /*!
   * Gives the position where the inlet can be found: 'top' or 'bottom'
   */
  string m_inlet_pos_str;

//  /*!
//   * The length of the inlet
  //   */
//  double m_inlet_length;

//  /*!
//   * The normal vector of the inlet
  //   */
//  point_t m_inlet_normal;
  
  /*!
   * Scale factor for the boundary
   */
  point_t m_scale;
  
  /*!
   * Is the STL ascii or binary? Has to be specified by the user.
   */
  string m_format;

//  /*!
//   * The base surface of the inlet is stored in this WallContainer
//   * so that ParticleCreatorInlet can determine the volume of 
//   * the inlet correctly.
  //   */
//  WallContainer m_inlet_base_surface;

  /*!
   * Initialize the property list
   */
  void init();
	
public:
  /*!
   * Constructor
   * @param phase Pointer to the parent \a Phase object this \a Boundary belongs to
   */
  BoundarySTL(Phase *phase);

  /*!
   * Destructor
   */
  virtual ~BoundarySTL();
	
  /*!
   * Initialize cell subdivision
   */
  virtual void setup(Simulation* sim, ManagerCell *mgr);

  /*!
   * Read and scale the boundary
   */
  virtual void setup();

//  /*!
//   * Return the length of the inlet
  //   */
//  double inletLength() const {
//    return m_inlet_length;
//  }
  
//  /*!
//   * Return the normal vector of the inlet
  //   */
//  const point_t &inletNormal() const {
//    return m_inlet_normal;
//  }

//  /*!
//   * Return the base surface of the inlet
  //   */
//  WallContainer &inletBaseSurface() {
//    return m_inlet_base_surface;
//  }
};

#endif
