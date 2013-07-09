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



#ifndef __F_DPDE_H
#define __F_DPDE_H 

#include "f_with_rng.h"
#include "pairdist.h"
#include "simulation.h"
#include "manager_cell.h"
#include "function_fixed.h"
#include "integrator_energy.h"
#include "weighting_function.h"


using namespace std;


/*!
 * Implementation of DPDE according to
 * J. B. Avalos and A. D. Mackie, Europhys. Lett. 40, 141-146 (1997) and
 * P. Espanol, Europhys. Lett. 40, 631-636 (1997)
 */
class FDPDE : public FWithRng
{
protected:
  /*!
   * For passing information to the threads.
   */
  size_t m_force_index;

  /*!
   * The dissipation constant.
   */
  double m_dissipation;

  /*!
   * The noise amplitude.
   */
  double m_noise;

  /*!
   * Reciprocal of the sqrt of the time step. For faster calculation.
   */
  double m_sqrt_r_dt;

  /*!
   * Integrators
   */
  pair<IntegratorEnergy*, IntegratorEnergy*> m_ie;

  /*!
   * Stores the index of the force for the force factors.
   */
  pair<size_t, size_t> m_eforce_offset;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   */
  FDPDE(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~FDPDE();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  /*!
   * Computes the force on the two particles.
   * @param force_index The index within the force history
   */
  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

  /*!
   * Setup shortly before simulation starts
   */
  virtual void setup();
};

#endif
