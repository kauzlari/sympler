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



#ifndef __F_PAIR_VELS_H
#define __F_PAIR_VELS_H 

#include "f_pair_arbitrary.h"
#include "pairdist.h"
#include "manager_cell.h"
#include "pair_creator.h"
#include "function_pair.h"
#include "weighting_function.h"


using namespace std;

/*!
 * This is a completely general pair force FPV on a velocity-vector v_i such that:
 * dv_i = FPV*dt
 *      = Sum_j(particleFactor(i)_ij*pairFactor_ij)*dt
 * where pairFactor_ij includes the (anti-)symmetric contributions of the pair ij,
 * and particleFactor(i)_ij the non-symmetric ones for particle i.
 * Note that the expression has to give a vector in the end,
 * otherwise the compiler will report an error.
 */
class FPairVels : public FPairArbitrary
{
protected:

/*!
 * Initialise the property list
 */
  void init();

public:
 /*!
 * Constructor
 * @param simulation The \a Simulation object the force belongs to
 */
  FPairVels(Simulation *simulation);

  /*!
  * Destructor
  */
  virtual ~FPairVels();

  /*!
   * Setup this force, mainly the slots in memory
   */
  virtual void setup();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  /*!
   * Compute the force
   * @param force_index The index for the memory slot to save the current force
     */
  virtual void computeForces(int force_index);

  /*!
   * Compute the force
   * @param force_index The index for the memory slot to save the current force
   * @param pair is the current Pairdist, the force acts on.
   */
#ifndef _OPENMP
void computeForces(Pairdist* pair, int force_index)
#else
void computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{

  /*fixme!!! there are weighting functions which use the particle position as
    argument but they are not yet usable here; that's why the following dirty hack*/
  if (this->m_cutoff > pair->abs()) 
    {                                    
      point_t temp;

      this->m_pairFactor(&temp, &(*pair));
            
      point_t fi;
      point_t fj;
      
      // compute the particle-expressions
      this->m_1stparticleFactor(&fi, &(*pair));
      this->m_2ndparticleFactor(&fj, &(*pair));
            
      // loop necessary because operator* of math_vector_t does scalar product
      for(size_t i = 0; i < SPACE_DIMS; ++i)
	{
	  fi[i] *= temp[i];
	  fj[i] *= temp[i];
	}
            
#ifndef _OPENMP
       if (pair->actsOnFirst()) {
	 pair->firstPart()->force[force_index] += fi;
       }
	 
#else

   if (pair->actsOnFirst()) {
     for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
	      (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi[_i];
     }
   }
#endif
              
#ifndef _OPENMP
	 if (pair->actsOnSecond()) {
	   pair->secondPart()->force[force_index] += m_symmetry*fj;
	 }
#else
     if (pair->actsOnSecond()) {
	     for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
	       (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += m_symmetry*fj[_i];
	     }
     }
#endif 
       
    }                                                                      
}



  /*!
   * Will throw an exception
   * @param force_index The index for the memory slot to save the current force
   * @param part is the current Particle this force acts on.
   */
  virtual void computeForces(Particle* part, int force_index);

};

#endif
