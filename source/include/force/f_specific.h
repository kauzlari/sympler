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



#ifndef _F_SPECIFIC_H
#define _F_SPECIFIC_H

#include "f_particle.h"
#include "list"

using namespace std;

/* --- Fspecific --- */

class Simulation;

typedef list<Particle*> SpecificParticleList;
typedef list<Particle*>::iterator SpecificParticleListItr;

/*!
 * Implementation of a constant per particle force
 * acts on specific particles and not a species.
 */
class Fspecific : public GenF
{
  public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
	Fspecific(Simulation *simulation);

#ifdef _OPENMP
       virtual void setForceSlots(Integrator* intr, int thread_no) {}
#endif

  /*!
   * Destructor
   */
	virtual ~Fspecific();
	virtual void computeForces(int force_index);

#ifndef _OPENMP
        virtual void computeForces(Pairdist* pair, int force_index);
#else
        virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//         virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

        virtual void computeForces(Particle* part, int force_index);

 protected:
  /*!
   * The force vector
   */
  point_t m_force;

  /*!
   */
  string m_species;

  /*!
   * The colour index that corresponds to the species \a m_species
   */
  size_t m_colour;

  /*!
   * Inititialize the property list.
   */
  void init();
  /*!
   * The particles the force acts on
   */
  SpecificParticleList m_SpecificParticleList;
  virtual void setup();
  public:
  /*!
  * Adds a particle pair to the connected list
   */
  void addParticleToForce(Particle *p);

  };
#endif
