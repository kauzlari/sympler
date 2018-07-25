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



#ifndef __F_GENERIC_H
#define __F_GENERIC_H 

#include "f_pair_wf.h"
#include "pairdist.h"
#include "function_fixed.h"
#include "simulation.h"
#include "manager_cell.h"
#include "integrator_energy.h"
#include "weighting_function.h"
#include "pca_energy_entropy.h"


using namespace std;


/*!
 * Implementation of a fluid particle model with degrees of freedom
 * r - position, v - velocity, V - volume, and e - internal energy
 */
class FGeneric : public FPairWF
{
protected:
  /*!
   * Stores the index of the force for the force factors.
   */
  pair<size_t, size_t> m_eforce_offset[FORCE_HIST_SIZE];

  /*!
   * Stores the index of the temperature
   */
  pair<size_t, size_t> m_temperature_offset;

//  /*!
//   * For passing information to the threads.
  //   */
//   size_t m_force_index;

  /*!
   * Cutoff of the density
   */
  double m_density_cutoff;

  /*!
   * The derivative of the excess local energy with respect to the local density
   */
  FunctionFixed m_de_dn;

  /*!
   * Pressure = Derivative of entropy with respect to n * T
   */
  FunctionFixed m_Tds_dn;

  /*!
   * Pointer to the cache calculating de_dn and Tds_dn and storing the information
   * for each particle
   */
  ParticleCacheEnergyEntropy *m_pcee;

  /*!
   * Surface entropy for defining contact angles
   */
  double m_surface_entropy;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   */
  FGeneric(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~FGeneric();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  /*!
   * The actual force computation
   * @param force_index The slot in the force history table
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
   * Register the particle cache, etc...
   */
  virtual void setup();
};

#endif
