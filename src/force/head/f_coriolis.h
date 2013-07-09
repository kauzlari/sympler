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



#ifndef __F_CORIOLIS_H
#define __F_CORIOLIS_H

#include "f_particle.h"

#include "function_fixed.h"

using namespace std;

/* --- FCoriolis --- */

class Simulation;

/*!
 * Coriolis force on reach particle, i.e. f = w x v
 */
class FCoriolis : public FParticle
{
protected:
  /*!
   * The rotation axis normalized to unit length
   */
  point_t m_rotation_axis;

  /*!
   * The rotation frequency
   */
  FunctionFixed m_frequency;

  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  FCoriolis(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~FCoriolis();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no) {}
#endif

  virtual void setup();

  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif
};

#endif
