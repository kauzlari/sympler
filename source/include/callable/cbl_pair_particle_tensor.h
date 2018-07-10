/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#ifndef __CBL_PAIR_PARTICLE_TENSOR_H
#define __CBL_PAIR_PARTICLE_TENSOR_H

#include "pair_creator.h"
#include "threads.h"

#include "cbl_pair_particle_arbitrary.h"

/*!
 * Callable computing completely user-defined tensor properties for 
 * each particles based on pair summation
 */
class CblPairParticleTensor : public CblPairParticleArbitrary
{
  
 protected:

  /*!
   * Initialise the property list
   */
  virtual void init();

  
 public:

  /*!
   * Constructor for the \a Node hierarchy
   */
  CblPairParticleTensor(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~CblPairParticleTensor();

  /*!
   * Compute the user defined expression for pair \a pD
   * FIXME: Think how to test or make testable
   * @param pD \a Pairdist whose contribution we calculate
   */
  virtual void call(size_t timestep) {

    Phase* phase = ((Simulation*) m_parent) -> phase();
    
    // will do nothing if already up to date
    phase->pairCreator()->createDistances();

    tensor_t temp, tempFirst, tempSecond;
    
    FOR_EACH_PAIR__PARALLEL
      (CblPairParticleTensor,
       m_cp,
       if(pair->abs() < self->m_cutoff) {

	 // compute the pair-expression
	 m_function(&temp, pair);
	 
	 Particle* first = pair->firstPart();
	 Particle* second = pair->secondPart();
	 
	 if(pair->actsOnFirst()) {
	   
	   m_1stparticleFactor(&tempFirst, pair);
	   
	   for(size_t i = 0; i < SPACE_DIMS; ++i)
	     for(size_t j = 0; j < SPACE_DIMS; ++j)
	       tempFirst(i, j) *= temp(i, j);
	   
	   first->tag.tensorByOffset(m_slots.first) += tempFirst;
	   
	 }
	 
	 if(pair->actsOnSecond()) {
	   
	   m_2ndparticleFactor(&tempSecond, pair);
	   
	   for(size_t i = 0; i < SPACE_DIMS; ++i)
	     for(size_t j = 0; j < SPACE_DIMS; ++j)
	       tempSecond(i, j) *= temp(i, j);
	   
	   second->tag.tensorByOffset(m_slots.second) += m_symmetry*(tempSecond);
	   
	 }
	 
       } // end of if(pair->abs() < self->m_cutoff)
       );
    
  }
  
  /*!
   * Setup this \a Callable
   */
  virtual void setup();

};

#endif
