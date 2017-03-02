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




#ifndef __INTEGRATOR_VELOCITY_VERLET_H
#define __INTEGRATOR_VELOCITY_VERLET_H

#include "integrator_position.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorVelocityVerlet ----

/*!
 * Modified Velocity-Verlet integrator for the positions and velocities
 * See: R. D. Groot and P. B. Warren, J. Chem. Phys. 107, 4423-4435 (1997)
 */

class IntegratorVelocityVerlet: public IntegratorPosition
{
protected:
  /*!
   * Initialize the property list
   */
  void init();
  /*!
   * Mass for the species this integrator works for.
   */
  double m_mass;
  /*!
   * Lambda of Warren and Groot paper
   */
  double m_lambda;
  /*!
   * Difference to lambda = 1/2 (the standard Velocity-Verlet)
   */
  double m_lambda_diff;



public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorVelocityVerlet(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorVelocityVerlet();

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
   * Integration of the position
   */
  virtual void integratePosition(Particle* p, Cell* cell);

  /*!
   * Prediction of the velocity
   */
  virtual void integrateVelocity(Particle* p);

  /*!
   * Solves the equation that checks for hits
   */
  virtual void solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t
&force, vector<double>* results);

  /*!
   * Checks which of the times (in the time vector) is the actual hit position. The function will
be used in WallTriangle
   */
  virtual void hitPos(/*WallTriangle* wallTriangle, */double dt, const Particle* p, point_t &hit_pos,
const point_t &force);

#ifdef _OPENMP
  virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

#endif

};


#endif


