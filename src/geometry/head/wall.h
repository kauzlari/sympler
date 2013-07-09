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



#ifndef __WALL_H
#define __WALL_H

#include <string>
#include <ostream>

#include "reflector.h"
#include "smart_enum.h"
#include "node_one_child.h"


class WallContainer;

//---- class Wall -----

/*!
 * Base class of a wall piece.
 */
class Wall: public NodeOneChild
{
 protected:
  /*!
   * Information for export to VTK: Cell type
   */
  int m_vtk_cell_type;

  /*!
   * Information for export to VTK: Number of vertices
   */
  int m_vtk_n_vertices;

  /*!
  * stores the periodicity of the \a Boundary where this \a Wall belongs to
  */
  bool_point_t m_periodicity;

  /*!
  * stores the box size of the \a Boundary where this \a Wall belongs to
  */
  point_t m_boxSize;
    
  /*!
   * The reflector assigned to this wall
   */
  Reflector *m_reflector;

  /*!
   * Does nothing
   */
  virtual Node *instantiateChild(const string &name);

 public:
  /*!
   * Constructor
   * @param container The \a WallContainer this wall belong to
   * @param reflector The \a Reflector for this wall
   */
  Wall(WallContainer *container, Reflector *reflector);

  /*!
   * Destructor
   */
  virtual ~Wall();
	
  /*!
   * Setup the \a Wall
   */
  virtual void setup();
  
  /*!
  * Set \a m_periodicity \a m_boxSize and 
  */
  virtual void setBoundaryData(bool_point_t per, point_t size);
  
   /*!
   * Bounding box of this wall piece. 
   * Needed to find the total extend of the WallContainer.
   */
  virtual cuboid_t boundingBox() = 0;

  /*!
   * Notify the wall that the position of one or more vertices has changed.
   */
  virtual void vertexChanged() = 0;
	
  /*!
   * Intersection with a line
   * @param c Line to intersect with
   * @param hit_pos Point of intersection
   */
  virtual bool intersects(const line_t &c, point_t &hit_pos) const = 0;

  /*!
   * Intersection with a line
   * @param from Starting point of the line to intersect with
   * @param dir direction of the line
   * @param dist resulting perpendicular distance of 'from' to the wall 
   */
  virtual bool intersects(const point_t& from, const point_t& dir, double& dist) const = 0;
      
  /*!
   * Intersection with a line which is parallel to one of the coordinate axises
   * @param l Line to intersect with
   * @param dir Coordinate axis
   * @param hit_pos Point of intersection
   */
  virtual bool intersectsParallelToDir(const line_t &l, int dir, point_t &hit_pos) const = 0;

  /*!
   * Intersection with a cuboid
   * @param c Cuboid to intersect with
   */
  virtual bool intersects(const cuboid_t &c) const = 0;

  /*!
   * Check whether a particle will hit this wall during this time step
   * @param p Check for this particle
   * @param force Force on this particle
   * @param t_traveled The time the particle travelled until the hit occured
   * @param hit_pos The position where the particle hit the wall
   */
  virtual bool hit(const Particle *p, const point_t &force, double &t_traveled, point_t &hit_pos, class IntegratorPosition* integratorP) = 0;

  /*!
   * distance to plane, the wall is in
   */
  virtual double distToPlane(const point_t& p) const = 0;
  
  /*!
   * Is this wall close than \a range to position \a p?
   * @param range Distance from this wall
   * @param p Position to check distance from
   */
  virtual bool isInRange(const double& range, const point_t& p) const = 0;
			
  /*!
   * Returns a normal vector. Fixme!!! This should be position depedent generally.
   */
  virtual const point_t& normal() const = 0;

  /*!
   * Returns an in plane vector. Fixme!!! This should be position depedent generally.
   */
  virtual const point_t& inPlane() const = 0;
	
  /*!
   * Returns the reflector.
   */
  virtual Reflector *reflector() {
    return m_reflector;
  }
	
  /*!
   * Return a string identifier for this wall
   */
  virtual string toString() const = 0;
	
  /*!
   * distance of a point to the wall (however it may be defined for concrete walls)
   * @param p Distance to this point
   */
  virtual double distanceTo(const point_t& p) const = 0;

  
  /*!
   * Check, whether the two walls are the same
   */
  virtual bool operator==(const Wall& wall) const = 0;

  /*!
   * Returns the VTK cell type
   */
  int vtkCellType() {
    return m_vtk_cell_type;
  }

  /*!
   * Returns the number of vertices needed for this VTK cell type.
   */
  int vtkNVertices() {
    return m_vtk_n_vertices;
  }
  
  /*!
   * Write this wall to VTK
   * @param s Stream to write to
   */
  virtual void toVTK(ostream &s) {
  }
};



//---- Factories ----

class Wall_Factory: public SmartEnum<Wall_Factory>
{
public:
    virtual Wall *instantiate(WallContainer *wall_container, Reflector *reflector) const = 0;
    
protected:
    Wall_Factory(const string &name)
    : SmartEnum<Wall_Factory>(name) { }
};


template <class T>
class Wall_Register: public Wall_Factory
{
public:
    Wall_Register(const string &name)
    : Wall_Factory(name) { }
    
    virtual Wall *instantiate(WallContainer *wall_container, Reflector *reflector) const {
        return new T(wall_container, reflector);
    }
};

#endif
