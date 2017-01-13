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


//
// C++ Implementation: integrator_rigidbodyMC
//
// Description: This integrator performs random moves of a rigid body for all its
// six degrees of freedom
//


#ifndef __INTEGRATOR_RIGIDBODYMC_H
#define __INTEGRATOR_RIGIDBODYMC_H

#include "integrator.h"
#include "random.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorRigidBodyMC----

/*!
 * C++ Implementation: integrator_rigidbodyMC
 *
 * Description: This integrator performs random moves of a rigid body for all its
 * six degrees of freedom
 */

class IntegratorRigidBodyMC: public Integrator
{
protected:
  /*!
   * The name of the integrator displacement
   */
  string m_displacement_name;

  /*!
   * The symbol (short name) of the integrator displacement
   */
  string m_displacement_symbol;
  /*!
   * The tag offset of the integrator displacement 
   */
  size_t m_displacement_offset;

  /*!
   * The tag offset of the last move
   */
  size_t m_lastmove_offset;
  /*!
   * Declare the random number generator
   */
    RandomNumberGenerator m_rng;
  /*!
   * Initialize the property list
   */
  void init();
  /*!
   * first time rotation matrix calculation needed, this boolean must be set to true
   */
  bool m_first_sweep;
  /*!
   *  normalize quaternions after timestep?
   */
  bool m_normalize;
  /*!
   * reject the move?
   */
  bool m_reject;
  /*!
   * the maximum translational and rotational displacements
   */
  double m_DeltaR;
  double m_PHIMax;
  /*!
   * the temperature we start the MC with
   */
  double m_kT;
  /*!
   * Intervals for writing the acceptance rate
   */
  int m_wint;
  /*!
   * The name of total energy.
   */
  string m_Etot_name;
  /*!
   * The tag offset of the total energy
   */
  size_t m_Etot_o;
  /*!
   * The doubles for the new and the old total energy.
   */
  double EtotNew;
  double EtotOld;
  /*!
   * The names of the quaternion components.
   */
  string m_q0_name;
  string m_q1_name;
  string m_q2_name;
  string m_q3_name;
  /*!
   * The symbols (short names) of the quaternion components
   */
  string m_q0_symbol;
  string m_q1_symbol;
  string m_q2_symbol;
  string m_q3_symbol;
  /*!
   * The tag offset of the quaternion components
   */
  size_t m_q0_offset;
  size_t m_q1_offset;
  size_t m_q2_offset;
  size_t m_q3_offset;
  /*!
   * The tag offset of old quaternion components
   */
  size_t m_q0Old_o;
  size_t m_q1Old_o;
  size_t m_q2Old_o;
  size_t m_q3Old_o;
    /*!
   * The name of the rotation matrix
   */
  string m_rotmat_name;
  /*! 
   * The symbol (short name) of the rotation matrix
   */
  string m_rotmat_symbol;
  /*!
   * The tag offset of the rotation matrix
   */
  size_t m_rotmat_offset;
  /*!
   * The name of the transpose of the rotation matrix
   */
  string m_rotmatT_name;
  /*! 
   * The symbol (short name) of the transpose of the rotation matrix
   */
  string m_rotmatT_symbol;
  /*!
   * The tag offset of the transpose of the rotation matrix
   */
  size_t m_rotmatT_offset;
  /*!
   * Time step
   */
  double m_dt;

  /*!
   * The maximum velocity for each MC move
   */
  double m_velmax;
  /*!
   * The number of moves m_moves and the number of accepted moves
   * m_accept
   */
  size_t m_moves, m_accept;
public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorRigidBodyMC(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorRigidBodyMC();

  /*!
    *    * Register the field and the force of the field with the \a Particle
    */ 
    virtual void setup();

  /*!
    *    * Calculate the rotation matrix
    */ 
    virtual void setupAfterParticleCreation();

   /*
      Called right before the simulation will start
      virtual void isAboutToStart();
   */

  /*!
   * Calculation of the new orientation
   */
  virtual void integrateStep1();

  /*!
   * If used forces are saved in a tag, protect and unprotect them as needed
   */
    virtual void unprotect(size_t index){
    }

  /*!
   * Calculation of the new angular velocity
   */
  virtual void integrateStep2();

  /*!
   * Integration of the position
   */
  virtual void integratePosition(Particle* p, Cell* cell);

  /*!
   * Prediction of the velocity
   */
  virtual void integrateVelocity(Particle* p){
  }
  /*!
   * Solves the equation that checks for hits
   */
  virtual void solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, 
      					const point_t &force, vector<double>* results);

  /*!
   * Checks which of the times (in the time vector) is the actual hit position. 
   * The function will be used in WallTriangle
   */
  virtual void hitPos(/*WallTriangle* wallTriangle, */double dt, const Particle* p, 
      					point_t &hit_pos, const point_t &force);

#ifdef _OPENMP
virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

  /*!
   * Returns the number of doubles this Integrator saves in a particle tag.
   */
    virtual int numCopyDoubles() {
      return SPACE_DIMS;
    }


#endif

};


#endif


