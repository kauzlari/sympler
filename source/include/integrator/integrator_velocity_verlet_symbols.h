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




#ifndef __INTEGRATOR_VELOCITY_VERLET_SYMBOLS_H
#define __INTEGRATOR_VELOCITY_VERLET_SYMBOLS_H

#include "integrator.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorVelocityVerletSymbols ----

/*!
 * Modified Velocity-Verlet integrator for the positions and velocities
 * See: R. D. Groot and P. B. Warren, J. Cham. Phys. 107, 4423-4435 (1997). 
 * This \a Intgrator integrates positions and velocities stores in the \a Particle 's tag, keeping the original positions and velocities untouched. Hence, this class may not inherit from \a IntegratorPosition
 */

class IntegratorVelocityVerletSymbols: public Integrator
{
protected:

  /*!
   * Full name of the velocity, usable as attribute in other modules
   */
  string m_vel_name;

  /*!
   * Full name of the position, usable as attribute in other modules
   */
  string m_pos_name;

  /*!
   * Symbol name for the velocity, usable in algebraic expressions
   */
  string m_vel_symbol;

  /*!
   * Symbol name for the position, usable in algebraic expressions
   */
  string m_pos_symbol;

  /*!
   * Memory offset in the \a Particle tag for the velocity, usable in algebraic expressions
   */
  size_t m_vel_offset;

  /*!
   * Memory offset in the \a Particle tag for the position, usable in algebraic expressions
   */
  size_t m_pos_offset;

  /*!
   * The tag offsets of the force on this tensor field.
   */
    size_t m_force_offset[FORCE_HIST_SIZE];

  /*!
   * The indices of the force attribute of the tag of the force on this tensor field.    
   */
    size_t m_fAttr_index[FORCE_HIST_SIZE];

  /*!
   * Lambda of Warren and Groot paper
   */
  double m_lambda;

  /*!
   * Difference to lambda = 1/2 (the standard Velocity-Verlet)
   */
  double m_lambda_diff;

  /*!
   * Timestep
   */
  double m_dt;

  /*!
   * Timestep divided by particle mass
   */
  double m_dt_div_mass;

  /*!
   * Half of timestep divided by particle mass
   */
  double m_dt_div2_mass;


  /*!
   * Initialize the property list
   */
  void init();
  /*!
   * Mass for the species this integrator works for.
   */
  double m_mass;


public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorVelocityVerletSymbols(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorVelocityVerletSymbols();



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

  /*!
   * Setup this \a Integrator
   */
  virtual void setup();

  /*!
   * If used forces are saved in a tag, protect and unprotect them as needed
   */
  virtual void unprotect(size_t index);

#ifdef _OPENMP
  virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

  /*!
   * Returns the number of double-copies this Integrator saves in a particle tag.
   */
  virtual int numCopyDoubles() {
    return SPACE_DIMS;
  }

#endif

};


#endif


