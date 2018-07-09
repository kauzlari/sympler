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



#ifndef __REFLECTOR_H
#define __REFLECTOR_H

#include "node.h"
#include "general.h"
#include "particle.h"
#include "smart_enum.h"
#include "random.h"

// class Wall;
class Boundary;

/*!
 * \a Reflector s define what happens if a particle hits a hard wall
 */
class Reflector: public Node
{
protected:
  /*!
   * String identifier for this reflector
   */
  string m_reflector_name;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param wall The wall this constructor belongs to. Fixme!!! Obsolete?
   */
  Reflector(/*Wall *wall*/ Boundary* boundary);

  /*!
   * Destructor
   */
  virtual ~Reflector();
  
  /*!
   * Return the name of this reflector
   */
  string name() {
    return m_reflector_name;
  }
    
  /*!
   * Reflect a particle
   * @param p Particle to reflect
   * @param hit_pos Position, where particle hit the wall
   * @param normal Normal vector to the wall
   * @param in_plane A vector perpendicular to the normal vector
   */
  virtual void reflect(Particle *p, point_t& r, point_t& v, const point_t &hit_pos, const point_t &normal, const point_t &in_plane) = 0;

  
  /*!
   * A counter for the total number of hits on the walls.
   */
  static size_t s_n_hits;
};



//---- Factories ----

class Reflector_Factory: public SmartEnum<Reflector_Factory>
{
public:
   virtual Reflector *instantiate(/*Wall *wall*/ Boundary* boundary) const = 0;
    
protected:
    Reflector_Factory(const string &name)
    : SmartEnum<Reflector_Factory>(name) { }
};


template <class T>
class Reflector_Register: public Reflector_Factory
{
public:
    Reflector_Register(const string &name)
    : Reflector_Factory(name) { }
    
    virtual Reflector *instantiate(/*Wall *wall*/ Boundary* boundary) const {
        return new T(/*wall*/ boundary);
    }
};

#endif
