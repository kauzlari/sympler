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



#ifndef __F_WALL_REPULSION_H
#define __F_WALL_REPULSION_H 

#include "f_particle.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


//---- FWallRepulsion ----

/*!
 * Simple repulsive wall force which drops linearly from \a m_force at the wall
 * to zero at the cut-off distance.
 */
class FWallRepulsion : public FParticle
{
protected:
  /*!
   * The cut-off distance from the wall
   */
  double m_cutoff;

  /*!
   * The inverse of the cut-off distance from the wall
   */
  double m_rcinv;

  /*!
   * Direction in where to find the walls
   */
  int m_wall_dir;

  /*!
   * Position of the wall to the left in the wall direction
   */
  double m_left_wall;

  /*!
   * Position of the wall to thr right in the wall direction
   */
  double m_right_wall;

  /*!
   * Magnitude of the force directly at the wall
   */
  double m_force;

  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  FWallRepulsion(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~FWallRepulsion();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no) {}
#endif

  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

  virtual void setup();
};

#endif
