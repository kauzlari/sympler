/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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



#ifndef __LJ_H
#define __LJ_H

#include "f_pair.h"
/* #include "gen_f_solo.h" */

using namespace std;

class LJ : public FPair

{
public:
  LJ(Simulation *simulation);
  virtual ~LJ() {}

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index)
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
  {
    // the inverse of the pair-distance
    double r = pair->abs();
    if (r < /*this->*/m_cutoff) {
      double r2i = 1/(r*r);
      double r6i = r2i*r2i*r2i;
      point_t f = 
	48*m_epsilon*(m_sigma_pow_12*r6i-m_half_sigma_pow_6)*r6i*r2i /*- m_shift_force*/ 
	* pair->cartesian();
      
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
      // m_c1, m_c2 have been set in setup()
      
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
      
    }
    
  }

#ifdef _OPENMP
  void setForceSlots(Integrator* intr, int thread_no);
#endif
		
  /*   virtual void computeEnergy(double& energy, group_t* for_groups); */

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
		
 protected:
  /*!
   * Well depth
   */
  double m_epsilon; 
  /*!
   * zero energy distance
   */
  double m_sigma; 
  /*!
   * helper
   */
  double m_half_sigma_pow_6;
  /*!
   * helper
   */
  double m_sigma_pow_12;

  /*!
   * First colour according to the corresponding \a ColourPair
   */
  size_t m_c1;

  /*!
   * Second colour according to the corresponding \a ColourPair
   */
  size_t m_c2;

  // shifting the force is currently not implemented
  /*double m_shift_force;*/ 
  /*   bool m_shift_potential; */
  /*   double m_shift_energy; */
  
  /*!
   * Initialise the \a PropertyList
   */
  void init();
};

#endif
