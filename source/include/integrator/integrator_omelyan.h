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
// C++ Implementation: integrator_omelyan
//
// Description: This integrator is based on the paper by Igor P. Omelyan,
// On the numerical integration of rigid polyatomics: The modified quaternion approach,
// Computers in Physics, Vol. 12, No. 1, p. 97 (1998).
// Careful, I use the notation Q=(q_0,\vect{q}) here as in Allen/Tidesley, Computer Simulation 
// of Liquids (Clarendon, Oxford, 1987)
//


#ifndef __INTEGRATOR_OMELYAN_H
#define __INTEGRATOR_OMELYAN_H

#include "integrator.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorOmelyan ----

/*!
 * C++ Implementation: integrator_omelyan
 *
 * Description: This integrator is based on the paper by Igor P. Omelyan,
 * On the numerical integration of rigid polyatomics: The modified quaternion approach,
 * Computers in Physics, Vol. 12, No. 1, p. 97 (1998).
 * Careful, I use the notation Q=(q_0,\vect{q}) here as in Allen/Tidesley, Computer Simulation 
 * of Liquids (Clarendon, Oxford, 1987)
 */

class IntegratorOmelyan: public Integrator
{
protected:
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
   * components of the inertial tensor for the species this integrator works for.
   */
  double m_J1;
  double m_J2;
  double m_J3;
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
   * The names of the components of the angular velocity.
   */
  string m_omega_name;
  /*!
   * The symbols (short names) of the  components of the angular velocity
   */
  string m_omega_symbol;
  /*!
   * The tag offset of the angular velocity components
   */
  size_t m_omega_offset;
  /*!
   * The tag offset of the torque components saved in particle
   */
  size_t m_torque_offset;
  /*
   * The tag offsets of the force on the angular velocity components. Even though this integrator
   * needs only one, we let it use as much as the most sophisticated integrator. This is
   * currently (2006/01/11) FORCE_HIST_SIZE=2
  size_t m_ox_force_offset[FORCE_HIST_SIZE];
  size_t m_oy_force_offset[FORCE_HIST_SIZE];
  size_t m_oz_force_offset[FORCE_HIST_SIZE];
   */
  /*!
   * The torque delivered by the user
   */
  size_t m_omega_force_offset[FORCE_HIST_SIZE];
  size_t m_omega_fAttr_index[FORCE_HIST_SIZE];
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

public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorOmelyan(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorOmelyan();

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
    virtual void unprotect(size_t index);

  /*!
   * Calculation of the new angular velocity
   */
  virtual void integrateStep2();
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


