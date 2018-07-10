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



#ifndef __GEN_F_SOLO_H
#define __GEN_F_SOLO_H 

#include "f_with_rng.h"
#include "pairdist.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;


#define FORCE_FACTOR_STR "force_factor:"



/*!
 * This class describes forces who are linearly dependent on the
 *  distance vector between two particles and have a cutoff.
 */
class GenFSolo : public FWithRng
{
protected:
  /*!
   * Are the force factors old and need they be recalculated?
   * This is set to true by the \a invalidate() method.
   */
  bool m_force_factors_old;

#if 0
  /*!
   * Stores the index of the force for the force factors.
   */
  int m_offset;
#endif

  /*!
   * The offset of the inverse of the particle distance within
   * the pair's tag.
   */
  size_t m_compute_ri_offset;

  /*!
   * For passing information to the threads.
   */
  int m_force_index;

  /*!
   * The inverse of the cut-off radius.
   */
  double m_rcinv;

  void init();

  /*!
   * Compute the force factor for pair \a pair. The force is given by
   * force_factor * weight * eij.
   * @param pair The pair for which to calculate the force factor
   */
  virtual double computeForceFactor(Pairdist *pair) = 0;

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  GenFSolo(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~GenFSolo();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index);
#endif

  virtual void setup();

  inline virtual void invalidate() {
    m_force_factors_old = true;
  }
};

#endif
