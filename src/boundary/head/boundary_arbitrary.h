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



#ifndef __BOUNDARY_ARBITRARY_H
#define __BOUNDARY_ARBITRARY_H

#include "boundary.h"
#include "wall_container.h"

//---- class BoundaryArbitrary -----

/*!
 * A boundary that is given by a triangulated surface.
 */
class BoundaryArbitrary : public Boundary
{
protected:	
  /*!
   * The name of the VTK to be written holding the geometry.
   */
  string m_geometry_filename;

  /*!
   * The container holding the triangulated surface
   */
  WallContainer *m_container;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param phase Pointer to the \a Phase object this \a Boundary belongs to
   */
  BoundaryArbitrary(Phase *phase);

  /*!
   * Destructor
   */
  virtual ~BoundaryArbitrary();
	
  /*!
   * Ask David
   */
  virtual const bool_point_t& periodicityFront() const {
    return m_container -> periodicityFront();
  }
    
  /*!
   * Ask David
   */
  virtual const bool_point_t& periodicityBack() const {
    return m_container -> periodicityBack();
  }

  /*!
   * Returns the volume of the bounding box of the geometry
   */
  virtual const double& cuboidVolume() const {
    return m_container -> cuboidVolume();
  }
		
  /*!
   * Returns the bounding box of the geometry
   */
  virtual const cuboid_t &boundingBox() const {
    return m_container->boundingBox();
  }

  /*!
   * Check whether position \a point is inside the boundary.
   * @param point Position to check if inside or outside
   */
  virtual bool isInside(point_t point);

  /*!
   * Check whether cuboid \a cuboid is partially inside the boundary.
   * @param cuboid Cuboid to check if inside or outside
   * @param range Fixme!!! Look up.
   */
  virtual bool isInside(cuboid_t cuboid, const double& range);

  /*!
   * Return whether there is a wall closer than \a range to position \a p.
   * @param range Wall distance
   * @param point Point from which to check for walls within distance \a range
   */
  virtual bool isInWallRange(const double& range, const point_t& point);
	
  /*!
   * Just calls Boundary::setup()
   */
  virtual void setup();

  /*!
   * Adjust the size if frozen particles are present
   */
  virtual void setup(Simulation* sim, ManagerCell *mgr);
};

#endif
