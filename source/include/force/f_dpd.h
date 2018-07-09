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



#ifndef __F_DPD_H
#define __F_DPD_H 

#include "f_with_rng.h"
#include "pairdist.h"
#include "simulation.h"
#include "manager_cell.h"
#include "function_fixed.h"
#include "weighting_function.h"


using namespace std;


//---- FDPD ----

/*!
 * Standard DPD force according to
 * P. Espanol and P. Warren, Europhys. Lett. 30, 191-196 (1995)
 */
class FDPD : public FWithRng
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
   * The temperature.
   */
  double m_temperature;

  /*!
   * reciprocal of the sqrt of the time step
   */
  double m_sqrt_r_dt;

  /*!
   * Initialization, i.e., register properties with the property list.
   */
  void init();

public:
  /*!
   * Constructor.
   */
  FDPD(Simulation *simulation);

  /*!
   * Destructor.
   */
  virtual ~FDPD();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  /*!
   * The actual force computation.
   */
  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index);
#endif

  /*!
   * \a setup() is called after \a init() has been called
   * for all objects.
   */
  virtual void setup();
};

#endif
