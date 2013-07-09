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



#ifndef __LJ_VC_H
#define __LJ_VC_H

#include "lennard_jones.h"

using namespace std;

class LJVC 
: public LJ
{
public:
  LJVC(Simulation *simulation);
  virtual ~LJVC() {}

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index)
#else
  virtual void computeForces(Pairdist* pair, int force_index, size_t thread_no)
#endif
  {
    if (pair->abs() < /*this->*/m_cutoff) {
      
      double r6i = pair->tag.doubleByOffset(m_compute_r6i_offset);
      double ri = pair->tag.doubleByOffset(m_compute_ri_offset);
      
      point_t f = 
	48*m_epsilon*(m_sigma_pow_12*r6i-m_half_sigma_pow_6)*r6i*ri*ri /*- m_shift_force*/ 
	* pair->cartesian();
      
#ifdef ENABLE_PTHREADS
      pair->firstPart()->lock();
      pair->secondPart()->lock();
#endif
      
      /* because of CONVENTION 2 in pairdist.h, actsOn*() will always return false for a
	 frozen particle and there is consequently no danger, e.g., of modifying a force
	 acting on a frozen particle, which does not exist in memory */
#ifndef _OPENMP     
      if (pair->actsOnFirst())
	pair->firstPart()->force[/* this->m_ */force_index] += f;
      if (pair->actsOnSecond())
	pair->secondPart()->force[/* this->m_ */force_index] -= f;
#else
      
      // FIXME!:
      // Are all these colour checks below really necessary?!?
      
      if (pair->actsOnFirst()) {
	for (int _i = 0; _i < SPACE_DIMS; ++_i) {
	  if (m_c1 == pair->firstPart()->c)
	    (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += f[_i];
	  else if (m_c2 == pair->firstPart()->c)
	    (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += f[_i];
	}
      }
      if (pair->actsOnSecond()) {
	for (int _i = 0; _i < SPACE_DIMS; ++_i) {
	  if (m_c1 == pair->secondPart()->c)
	    (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] -= f[_i];
	  else if (m_c2 == pair->secondPart()->c)
	    (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] -= f[_i];
	}
      }
      
#endif
      
#ifdef ENABLE_PTHREADS
      pair->secondPart()->unlock();
      pair->firstPart()->unlock();
#endif
    }
    
  }

  /*!
   * Setup this force
   */
  virtual void setup();

  /*!
   * Will throw an exception
   */
  virtual void computeForces(Particle* part, int force_index);

  /*!
   * Will throw an exception
   */
  virtual void computeForces(int force_index);
		
 private:
  /*!
   * Memeory offset to 1/r stored in the particle tag
   */
  size_t m_compute_ri_offset;
  
  /*!
   * Memeory offset to (1/r)^6 stored in the particle tag
   */
  size_t m_compute_r6i_offset;
    
  /*!
   * Initialise the \a PropertyList
   */
  void init();
};

#endif
