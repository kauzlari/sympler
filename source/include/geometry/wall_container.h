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



#ifndef __WALL_CONTAINER_H
#define __WALL_CONTAINER_H 

#include <list>
#include <vector>
#include <ostream>

#include "wall.h"
#include "reflector.h"
#include "vertex_list.h"


/*!
 * A \a WallContainer holds and manages a list of \a Wall s, which can right now be
 * only triangles.
 */
class WallContainer: public VertexList
{
 protected:
  /*!
   * List of the walls
   */
  list<Wall*> m_walls;

  /*!
   * Indicates whether the "lower end" of a direction is periodic or not
   */
  bool_point_t m_periodicFront;

  /*!
   * Indicates whether the "upper end" of a direction is periodic or not
   */
  bool_point_t m_periodicBack;

  /*!
   * The bounding box of the geometry
   */
  cuboid_t m_bounding_box;
	
  /*!
   * The volume occupied by the bounding box
   */
  double m_cuboidVolume;
	
  /*!
   * The default \a Reflector
   */
  Reflector *m_reflector;
	
 public:
  /*!
   * Constructor
   * @param reflector The default reflector
   */
  WallContainer(Reflector *reflector = NULL);

  /*!
   * Destructor
   */
  virtual ~WallContainer();

  /*!
   * Sets the periodicity of the "lower end" of all directions
   */
  virtual void setPeriodicityFront(const bool_point_t &periodic) {
    m_periodicFront = periodic;
  }

  /*!
   * Sets the periodicity of the "upper end" of all directions
   */
  virtual void setPeriodicityBack(const bool_point_t &periodic) {
    m_periodicBack = periodic;
  }
	
  /*!
   * Returns the periodicity of the "lower end" of all directions
   */
  virtual const bool_point_t& periodicityFront() const {
    return m_periodicFront;
  }	

  /*!
   * Returns the periodicity of the "upper end" of all directions
   */
  virtual const bool_point_t& periodicityBack() const {
    return m_periodicBack;
  }	
		
  /*!
   * Stretch the whole geometry by factor \a factor
   * @param factor Factor for stretch in x, y and z
   */
  virtual void stretchBy(const point_t &factor);

  /*!
   * Move all vertices by a certain distance
   * @param dist Distance to move vertices
   */
  virtual void moveVertices(const point_t& dist);
	
  /*!
   * Count the number of walls that intersect line \a l
   * @param l Line for intersection test
   */
  virtual int intersectionsGeneral(const line_t &l);

  /*!
   * Count the number of walls in \a walls that intersect line \a l, where
   * \a l is aligned parallel to axis \a dir
   * @param walls Walls for intersection test
   * @param l Line for intersection test
   * @param dir Direction the line is aligned to
   */
  static int intersections(list<Wall*> walls, const line_t &l, int dir);

  /*!
   * Check if point \a pos is inside the geometry
   * @param pos Check if this pos is inside
   */
  virtual bool isInside(const point_t &pos);

  /*!
   * Check if the cuboid \a cuboid is inside this geometry
   * @param cuboid Is this cuboid inside?
   * @param range Safety distance
   */
  virtual bool isInside(const cuboid_t &cuboid, const double& range);

  /*!
   * Is any wall closer than \a range to position \a p?
   * @param range Distance from this wall
   * @param point Position to check distance from
   */
  virtual bool isInWallRange(const double& range, const point_t& point);
		
  /*!
   * Add a new wall
   * @param wall Wall to add
   */
  virtual void addWall(Wall *wall) {
    m_walls.push_back(wall);
  }
	
  /*!
   * Return the bounding box of the geometry
   */
  virtual const cuboid_t &boundingBox() const {
    return m_bounding_box;
  }

  /*!
   * Notify all walls that at least one vertex has changed.
   */
  virtual void vertexChanged();
	
  /*!
   * Return all walls
   */
  virtual list<Wall*> &walls() {
    return m_walls;
  }

  /*!
   * Delete a wall
   * @param wall Delete the wall given by this iterator
   */
  virtual void deleteWall(list<Wall*>::iterator &wall) {
    m_walls.erase(wall);
  }
	
  /*!
   * Returns the default reflector
   */
  virtual Reflector *reflector() {
    return m_reflector;
  }
	
  /*!
   * Returns the volume within the bounding box
   */
  virtual const double& cuboidVolume() const {
    return m_cuboidVolume;
  }
	
  /*!
   * Update the bounding box information
   */
  virtual void updateBoundingBox();
	
  /*!
   * Export boundary VTK
   * @param filename Write to VTK with this filename
   */
  virtual void toVTK(string filename);

  /*!
   * Export boundary VTK
   * @param s Write to this stream
   */
  virtual void toVTK(ostream &s);
};

#endif
