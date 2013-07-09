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



#ifndef __F_PAIR_SCALAR_H
#define __F_PAIR_SCALAR_H 

#include "f_pair_arbitrary_wf.h"
#include "pairdist.h"
#include "manager_cell.h"
#include "function_pair.h"
#include "function_particle.h"
#include "weighting_function.h"

using namespace std;

/*!
 * Note: This is a completely general pair force FPS on a scalar s_i:
 *
 * ds_i = FPS*dt
 *      = particleFactor_i*Sum_j(pairFactor_ij*interpolant_ij)*dt
 *
 * where particleFactor_i is a vector and a sum of quantities related to particle i,
 *       pairFactor_ij is a vector and includes all pair contributions of the pair ij
 *       interpolant_ij represents the interpolation function
 *
 * note that the expression has to give a scalar in the end,
 * otherwise the compiler will report an error
 */
class FPairScalar : public FPairArbitraryWF
{
protected:
  // FIXME: Next two are general for each FPairArbitrary except for FPairVels. That's why it's still here  => refine hierarchy or remove FPairVels
  /*!
   * The name of the scalar this force acts on
   */
  string m_scalar_name;

  /*!
   * The offset of the forces on that scalar
   */
  pair<size_t, size_t> m_force_offset[FORCE_HIST_SIZE];

#if 0
  /*!
   * The pair contribution to the force
   */
  FunctionPair m_pairFactor;

  /*!
   * An additional factor for the 1st particle
   */
  FunctionParticle m_1stparticleFactor;

  /*!
   * An additional factor for the 2nd particle
   */
  FunctionParticle m_2ndparticleFactor;

  /*!
   * The expression for the particle factor
   */
  string m_particleFactorString;

  /*!
   * Symmetry of the force
   */
  int m_symmetry;
#endif

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   */
  FPairScalar(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~FPairScalar();

  /*!
   * Setup
   */
  virtual void setup();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  /*!
   * Compute the forces
   * @param force_index The slot in the force history to use
   */
  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif
};

#endif
