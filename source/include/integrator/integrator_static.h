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





#ifndef __INTEGRATOR_STATIC_H
#define __INTEGRATOR_STATIC_H

#include "integrator_position.h"
#include "function_fixed.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//---- IntegratorPosition ----


class IntegratorStatic: public IntegratorPosition
{
protected:

    /*!
   * User defined particle velocity with variables x, y, z, t: x-component
   */
  FunctionFixed m_velocity_x;

    /*!
   * User defined particle velocity with variables x, y, z, t: y-component
   */
  FunctionFixed m_velocity_y;

    /*!
   * User defined particle velocity with variables x, y, z, t: z-component
   */
  FunctionFixed m_velocity_z;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorStatic(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorStatic();

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
   * Integration of the velocity
   */
  virtual void integrateVelocity(Particle* p);

  /*!
   *  Setup the Integrator
   */
  virtual void setup();

#ifdef _OPENMP
  virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

#endif

  /*!
   * Solves the equation that checks for hits
   */
  virtual void solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t
&force, vector<double>* results);

  /*!
   * Checks which of the times (in the time vector) is the actual hit position. The function will
be used in WallTriangle
   */
  virtual void hitPos(const double& dt, const Particle* p, point_t &hit_pos,
const point_t &force);

};

#endif
