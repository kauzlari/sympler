/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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



#ifndef __INTEGRATOR_POSITION_H
#define __INTEGRATOR_POSITION_H

#include "integrator.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;
//class Wall;

//---- IntegratorPosition ----

/*!
 * Parent class for \a Integrator s integrating particle positions (and maybe other variables as well)
 * FIXME: We should make this an abstract class with pure virtual functions!
 */
class IntegratorPosition: public Integrator
{
protected:

  /*!
   * Timestep
   */
  double m_dt;

  /*!
   * Mass for the species this integrator works for.
   */
  double m_mass;

  /*!
   * Timestep divided by particle mass
   */
  double m_dt_div_mass;

  /*!
   * Half of timestep divided by particle mass
   */
  double m_dt_div2_mass;

//  /*!
//   * Index that points to the previous entry in the force history table
//   */
//   size_t m_other_force_index;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorPosition(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorPosition();

#ifdef _OPENMP
  /*!
   * Returns the number of double-copies this Integrator saves in a particle tag.
   */
  virtual int numCopyDoubles() {
    return SPACE_DIMS;
  }

  virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

#endif

  /*!
   * Initialize temporary fields and clear all forces
   */
  virtual void isAboutToStart();

  /*!
   * Estimation of the velocity
   */
  virtual void integrateStep1();

  /*!
   * Correction of the velocity
   */
  virtual void integrateStep2();

  /*!
   * Integrates the position of the particle according to the different Position Integrators.
   */
  virtual void integratePosition(Particle* p, Cell* cell);

  /*!
   * Integrates the velocities of the particle according to the different Position Integrators
   */
  virtual void integrateVelocity(Particle* p);

  /*!
   * Calculates the algorithm for checking the hit time that will be used in WallTriangle for hit-check
   */
  virtual void solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results);

  /*!
   * Chacks which of the times (in the time vector) is the actual hit position
   */
  virtual void hitPos(const double& dt, const Particle* p, point_t &hit_pos, const point_t &force);

  /*!
   * Return the time step
   */
  double dt() const {
    return m_dt;
  }

  /*!
   * Return timestep divided by particle mass
   */
  double dtDivMass() const {
    return m_dt_div_mass;
  }

  /*!
   * Return half of timestep divided by particle mass
   */
  double dtDiv2Mass() const {
    return m_dt_div2_mass;
  }
};

#endif
