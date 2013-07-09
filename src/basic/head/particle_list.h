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



#ifndef __PARTICLE_LIST_H
#define __PARTICLE_LIST_H

#include "particle.h"
#include "smart_pointer.h"

/*--- ParticleList ---*/

/*
  class vector_particle_sp: public SmartPointer< vector<Particle> >
  {
  public:
  inline vector_particle_sp(): SmartPointer< vector<Particle> >() { }

  inline Particle &operator[](int i) {
  return (*m_value)[i];
  }

  inline const Particle &operator[](int i) const {
  return (*m_value)[i];
  }    
  };
*/


/* Use this as a pointer to a particle. */
//typedef list<Particle>::iterator particle_p;

//class ParticleList: public list<Particle>

/*!
 * List of \a Particles 's
 */
class ParticleList: public SmartList<Particle>
{
 public:
  /*!
   * Constructor
   */
  ParticleList();

  /*!
   * Destructor
   */
  virtual ~ParticleList();

  /*!
   * Calculate the center of mass velocity
   */
  virtual point_t centerOfMassVelocity();

  /*!
   * Rescale all velocity in order to find a velocity
   * distribution with a variance that corresponds to temperature
   * \a temperature.
   * @param temperature Temperature to scale to
   */
  virtual void scaleVels(double temperature);
};

#endif
