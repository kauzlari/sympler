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



// Boundary only for 3D !!!
// this Boundary has periodic boundary conditions at all six walls of a cuboid 

#ifndef __BOUNDARY_H
#define __BOUNDARY_H

#include "general.h"
/*#include "reflector.h"*/
#include "smart_enum.h"
/*#include "particle_creator.h"*/
#include "node_many_children.h"


//---- Forward declarations ----

class Phase;
class Simulation;
class Reflector;
class ParticleCreator;
class ManagerCell;
    
//---- Classes ----

using namespace std;


/*!
 * Base class for the description of simulation boundaries.
 */
class Boundary: public NodeManyChildren
{
protected:
  /*!
   * The boundary proposes a size for the bounding box of the simulation domain. 
    * It may be changed by a \a ParticleCreator
   */
  point_t m_proposedSize;

  /*!
   * Extend the simulation box in negative x, y, z direction to allow for
   * the addition of frozen particles?
   */
  bool_point_t m_frontFrame;

  /*!
   * Extend the simulation box in positive x, y, z direction to allow for
   * the addition of frozen particles?
   */
  bool_point_t m_endFrame;

  /*!
   * List of particle creators
   */
  vector<ParticleCreator*> m_pcList;

  /*!
   * Resizing the boundary using cutoff (without "skin")
   */
  double m_thickness;

  /*!
   * List of reflectors
   */
  vector<Reflector*> m_reflectors;

  virtual Node *instantiateChild(const string &name);

  /*!
   * Find a reflector by its name
   * @param name Name of the reflector to find
   */
  Reflector *findReflector(const string &name);
	
public:
  /*!
   * Default constructor. Should not be used.
   */
  Boundary();

  /*!
   * Constructor.
   * @param phase Pointer to the phase object
   */
  Boundary(Phase *phase);

  /*!
   * Destructor
   */
  virtual ~Boundary();

  /*!
   * Checks if there are any \a ParticleCreator s and \a Reflector s defined,
   * otherwise boil out.
   */
  virtual void setup();
  
  /*!
   * Do we have a periodic face at the "lower" end of the three directions?
   */
  virtual const bool_point_t& periodicityFront() const = 0;

  /*!
   * Do we have a periodic face at the "upper" end of the three directions?
   */
  virtual const bool_point_t& periodicityBack() const = 0;

  /*!
   * Return whether there is a wall closer than \a range to position \a p.
   * @param range Wall distance
   * @param p Point from which to check for walls within distance \a range
   */
  virtual bool isInWallRange(const double& range, const point_t& p) = 0;

  /*!
   * Return the bounding box of this geometry
   */
  virtual const cuboid_t &boundingBox() const = 0;
		
  /*!
   * Return the volume occupied by the bounding box of this geometry
   */
  virtual const double& cuboidVolume() const = 0;

  /*!
   * Check whether position \a point is inside the boundary.
   * The default implementation for isInside just uses m_size as a reference.
   * @param point Position to check if inside or outside
   */
  virtual bool isInside(point_t point);

  /*!
   * Return the size to resize the boundary (without "skin").
   */
  virtual double thickness() {
    return m_thickness;
  }

  /*!
   * Check whether cuboid \a cuboid is partially inside the boundary.
   * Fixme!!! The default implementation for this is crap. 
   * @param cuboid Cuboid to check if inside or outside
   * @param range Fixme!!! Look up.
   */
  virtual bool isInside(cuboid_t cuboid, const double& range);
    
  /*!
   * Now the boundaries are responsible to initialize the cell subdivision when
   * ManagerCell is being used. After setup is called the size and volume property
   * of the boundary have to remain fixed.
   * @param sim Pointer to the simulation object
   * @param mgr Pointer to the cell subdivision manager
   */
  virtual void setup(Simulation* sim, ManagerCell *mgr);

  /*!
   * Create particle initially by calling the \a ParticleCreator s
   */
  virtual void createParticles();

  /*!
   * Create additional particles during run by calling the \a ParticleCreator s
   */
  virtual void createMoreParticles();
};



//---- Factories ----

class Boundary_Factory: public SmartEnum<Boundary_Factory>
{
public:
  virtual Boundary *instantiate(Phase *phase) const = 0;

protected:
  Boundary_Factory(const string &name)
      : SmartEnum<Boundary_Factory>(name) { }
};


template <class T>
class Boundary_Register: public Boundary_Factory
{
public:
    Boundary_Register(const string &name)
    : Boundary_Factory(name) { }

    virtual Boundary *instantiate(Phase *phase) const;
};



//---- Inline functions ----

template <class T>
inline Boundary *Boundary_Register<T>::instantiate(Phase *phase) const
{
    return new T(phase);
}


#endif
