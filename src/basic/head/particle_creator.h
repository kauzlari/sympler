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



#ifndef __PARTICLE_CREATOR_H 
#define __PARTICLE_CREATOR_H 

#include "node.h"
#include "calc.h"
#include "smart_enum.h"
#include "particle_list.h"
#include "random.h"


class ParticleCalculator;

// struct dof_info_t {
//   size_t offset;
//   ParticleCalculator *pc;
// //   string name;
// };


//---- Forward declarations ----

class Boundary;

//---- Classes ----

/*!
 * \a ParticleCreator 's are in charge of creating particles at the beginning of the 
 * simulation, and, sometimes, during the course of the simulation. This is, for example,
 * necessary if one wants an inlet within the simulation domain where one needs to
 * insert particles continuously.
 */
class ParticleCreator: public Node
{
 protected:
  /*!
   * The size of the bounding box of the boundary
   */
  point_t m_box_size;

  /*!
   * Create particles of this species
   */
  string m_species;

  /*!
   * Create particles of this color
   */
  size_t m_colour;

  /*!
   * Initialize the property list
   */
  void init();

	
 public:

  /*!
   * Variable which is in charge of informing other files whether we create frozen particles or not.
   */
  static bool s_createFrozenParts;

  /*!
   * Constructor. Cannot be used.
   */
  ParticleCreator();

  /*!
   * Constructor
   * @param boundary The \a Boundary this \a ParticleCreator belongs to
   */
  ParticleCreator(Boundary *boundary);

  /*!
   * Destructor
   */
  virtual ~ParticleCreator();

  /*!
   * This function adjusts the \a size, \a frameRCfront and \a frameRCend parameters
   * depending on whether, the boundary needs to be resized, or frozen particles
   * need to be added outside the boundary.
   * Fixme!!! Ugly, rethink!
   * @param size The boundary size
   * @param frameRCfront Increase the boundary size in front
   * @param frameRCend Increase the boundary size at the end
   */
  virtual void adjustBoxSize
    (point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend);

  /*!
   * Create the initial particles.
   */
  virtual void createParticles() = 0;

  /*!
   * Randomly create a particle somewhere
   */
  virtual void createParticle(bool assign_to_cell = true) {
    throw gError("ParticleCreator::createParticle", "Not implemented for: '" + className() + "'");
  }

  /*!
   * Randomly create \a n particles somewhere before the first timestep
   * @param n Number of particles to create
   */
  virtual void createMoreParticlesAtStart(size_t n);

  /*!
   * Create particles during the run if the derived ParticleCreator permits that
   * When this function is called, the PC should know itself, how many and what to create
   */
  virtual void createMoreParticles() {}

  /*!
   * Randomly delete a particle from the creation region
   */
  virtual void deleteParticle() {
    throw gError("ParticleCreator::deleteParticle", "Not implemented for: '" + className() + "'");
  }

  /*!
   * Move particles to the phase if created off-line
   */
  virtual void flushParticles();
		
  /*!
   * Write configuration to XML file.
   */
  virtual ostream &write(ostream &s, int shift = 0) = 0;

  /*!
   * Register the species (if not existent) and get the appropriate color
   */
  virtual void setup();
};



//---- Factories ----

class ParticleCreator_Factory: public SmartEnum<ParticleCreator_Factory>
{
public:
    virtual ParticleCreator *instantiate(Boundary *boundary) const = 0;

protected:
    ParticleCreator_Factory(const string &name)
        : SmartEnum<ParticleCreator_Factory>(name) { }
};


template <class T>
class ParticleCreator_Register: public ParticleCreator_Factory
{
public:
    ParticleCreator_Register(const string &name)
        : ParticleCreator_Factory(name) { }
	
    virtual ParticleCreator *instantiate(Boundary *boundary) const;
};



//---- Inline functions ----

template <class T>
inline ParticleCreator *ParticleCreator_Register<T>::instantiate(Boundary *boundary) const
{
    return new T(boundary);
}


#endif
