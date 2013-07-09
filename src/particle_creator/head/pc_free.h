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



#ifndef __PARTICLE_CREATOR_FREE_H
#define __PARTICLE_CREATOR_FREE_H

#include <vector>

#include "particle_creator.h"
#include "random.h"

//---- Forward declarations ----

class Boundary;


//---- Classes ----


/*! 
* General interface for \a ParticleCreator s creating free particles. 
*/
class ParticleCreatorFree: public ParticleCreator
{
protected:
  
  /*!
 * Kinetic temperature of the particles to be created
   */
  double m_temperature;
  
  /*!
   * Square root of the kinetic temperature of the particles to be created
   */
  double m_temperature_sqrt;
  
  /*!
  * Boundary information
  */
  point_t m_box_size;

  /*!
  * Temporary storage of particles
  */
  map<int, ParticleList> m_particles;
	
  /*!
  * Have the particles already been saved for a new run using a \a ParticleCreatorFile ?
  */
  static bool toFileDone;

  /*!
  * Apply user-defined transformations of the particle positions (p.r)
  */
  virtual void transformPos(Particle &p) = 0;
      
  /*!
   * Apply user-defined transformations of the particle velocities (p.v)
   */
  virtual void transformVel(Particle &p) = 0;
  
  /*!
  * Initialise the property list
  */
  void init();
	
public:
  ParticleCreatorFree();
  ParticleCreatorFree(Boundary *boundary);
  virtual ~ParticleCreatorFree();

  virtual void setup();
  virtual void flushParticles();
  virtual void flushParticles(Particle** first_p);
  virtual ostream &write(ostream &s, int shift = 0);
};



// //---- Factories ----
// 
// class ParticleCreatorFree_Factory: public SmartEnum<ParticleCreatorFree_Factory>
// {
// public:
//     virtual ParticleCreatorFree *instantiate(Boundary *boundary) const = 0;
// 
// protected:
//     ParticleCreatorFree_Factory(const string &name)
//         : SmartEnum<ParticleCreatorFree_Factory>(name) { }
// };
// 
// 
// template <class T>
// class ParticleCreatorFree_Register: public ParticleCreatorFree_Factory
// {
// public:
//     ParticleCreatorFree_Register(const string &name)
//         : ParticleCreatorFree_Factory(name) { }
// 	
//     virtual ParticleCreatorFree *instantiate(Boundary *boundary) const;
// };
// 
// 
// 
// //---- Inline functions ----
// 
// template <class T>
// inline ParticleCreatorFree *ParticleCreatorFree_Register<T>::instantiate(Boundary *boundary) const
// {
//     return new T(boundary);
// }
// 

#endif
