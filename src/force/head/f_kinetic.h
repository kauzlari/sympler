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



#ifndef __F_KINETIC_H
#define __F_KINETIC_H 

#include "f_pair_wf.h"
#include "pairdist.h"
#include "manager_cell.h"
#include "function_fixed.h"
#include "weighting_function.h"


using namespace std;


/*!
 * Implementation of transport equation
 */

class FKinetic : public FPairWF
{
protected:
  /*!
   * Relaxation constant
   */
  double m_lambda;

  /*!
   * Will the local density be computed over all ColourPairs (CP) or only over the CP
   * corresponding to the chosen species?
   */
  bool m_oneProp;

  /*!
   * Name of the velocity correlation degree of freedom
   */
  string m_vv_name;

  /*!
   * Name of the local density to be used
   */
  string m_rhoSymbol;

  /*!
   * Offset for access to the density for each species
   */
  pair<size_t, size_t> m_density_offset;

  /*!
   * Offset for access to the velocity correlation for each species
   */
  pair<size_t, size_t> m_vv_offset;

  /*!
   * For passing information to the threads.
   * The current index in the force history.
   */
  pair<size_t, size_t> m_vv_force_offset;

  /*!
   * For the execution of parallel loops. So each thread knows which force index
   * is the current one.
   */
  int m_force_index;

  void init();

public:
  /*!
   * Constructor.
   * @param simulation Pointer to the simulation object
   */
  FKinetic(Simulation *simulation);

  /*!
   * Destructor.
   */
  virtual ~FKinetic();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
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
