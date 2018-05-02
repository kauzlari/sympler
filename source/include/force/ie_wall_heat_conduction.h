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



#ifndef __IE_WALL_HEAT_CONDUCTION_H
#define __IE_WALL_HEAT_CONDUCTION_H 

#include "function.h"
#include "f_particle.h"
#include "simulation.h"
#include "manager_cell.h"
#include "integrator_energy.h"


using namespace std;


//---- IEWallHeatConduction ----

/*!
 * Heat conduction term between a wall and the particles
 */
class IEWallHeatConduction : public FParticle
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
   * Position of the wall to the right in the wall direction
   */
  double m_right_wall;

  /*!
   * = 1/sqrt(dt)
   */
  double m_r_sqrt_dt;

  /*!
   * Direction parallel to the wall and perpendicular to \a m_perp_dir2
   */
  int m_perp_dir1;

  /*!
   * Direction parallel to the wall and perpendicular to \a m_perp_dir1
   */
  int m_perp_dir2;

  /*!
   * Temperature of the left wall
   */
  double m_left_temperature;

  /*!
   * Temperaturee of the right wall
   */
  double m_right_temperature;

  /*!
   * Reciprocal of the temperature of the left wall
   */
  double m_lTinv;

  /*!
   * Reciprocal of the temperature of the right wall
   */
  double m_rTinv;

  /*!
   * If \a m_spot_size is bigger than zero, only a square spot of "radius"
   * \a m_spot_size centered at (0, 0) on both walls will be heated.
   */
  double m_spot_size;

  /*!
   * Dissipation constant
   */
  double m_kappa;

  /*!
   * Noise amplitude
   */
  double m_alpha;

  /*!
   * Tag offset of the force on the internal energy
   */
  size_t m_eforce_offset[FORCE_HIST_SIZE];

  /*!
   * Integrator for the internal energy. Needed for temperature calculation.
   */
  IntegratorEnergy* m_ie;

  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  IEWallHeatConduction(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~IEWallHeatConduction();

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
