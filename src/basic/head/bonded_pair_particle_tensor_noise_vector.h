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


#ifndef __BONDED_PAIR_PARTICLE_TENSOR_NOISE_VECTOR_H
#define __BONDED_PAIR_PARTICLE_TENSOR_NOISE_VECTOR_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
/* #include "val_calculator_arbitrary.h" */
#include "bonded_pair_particle_calc.h"
#include "colour_pair.h"
#include "function_pair_fixed.h"

/*!
 * Class doing summation over bonded pairs and computing a stochastic vector from a noise tensor with independent random elements obtained from a Gaussian distribution with unit variance. 
 */
class BondedPairParticleTensorNoiseVector : public BondedPairParticleCalc
{
 protected:

  RandomNumberGenerator m_rng;

  bool m_randomize;

  //  /*!
  //   * The noise amplitude
  //   */
  //  double m_noise;
  
  //  /*!
  //   * = m_noise/sqrt(dt)
  //   */
  //  double m_noise_and_time;


  /*!
   * pair contribution to the noise. It knows the symbol "{dW}" for the matrix of independent Wiener increments 
   */
    FunctionPairFixed m_pairFactor;

   /*!
     * The mathematical expression for \a m_pairFactor
    */
    string m_pairFactorStr;
    
      /*!
     * An additional factor for the 1st particle. It knows the symbol "{dW}" for the matrix of independent Wiener increments 
       */
    FunctionPairFixed m_1stparticleFactor;

    /*!
     * An additional factor for the 2nd particle. It knows the symbol "{dW}" for the matrix of independent Wiener increments 
     */
    FunctionPairFixed m_2ndparticleFactor;

   /*!
     * The mathematical expression for \a m_1stparticleFactor
    */
    string m_1stPExpression;
    
   /*!
     * The mathematical expression for \a m_2ndparticleFactor
    */
    string m_2ndPExpression;


  /*!
   * Pointer to \a Controller to access the time step size
   */
  Controller* m_controller;
  
  /*!
   * Initialise the property list
   */
  virtual void init();
  
  /*!
   * Helper function for polymorphic copying
   */
  virtual ValCalculator* copyMySelf()
  {
    return new BondedPairParticleTensorNoiseVector(*this);
  }
  
  public:
   
    /*!
   * Constructor for the \a Node hierarchy
     */
    BondedPairParticleTensorNoiseVector(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~BondedPairParticleTensorNoiseVector();

#ifdef _OPENMP
    /*!
     * Merge the copies of all threads together
     */
    virtual void mergeCopies(ColourPair* cp, int thread_no);
#endif

    /*!
     * Compute the user defined expression for pair \a pair
     * @param pair \a Pairdist whose contribution we calculate
     */
#ifndef _OPENMP
        virtual void compute(Pairdist* pair) 
#else
        virtual void compute(Pairdist* pair, int thread_no)
#endif
        {

	  double noise[SPACE_DIMS][SPACE_DIMS];

	  for(size_t i = 0; i < SPACE_DIMS; ++i) {
	    for(size_t j = 0; j < SPACE_DIMS; ++j) {
 	      noise[i][j] = m_rng.normal(1);
	    }
	  }

	  double* noisePointer = &(noise[0][0]);
 
	  point_t temp;           

	  m_pairFactor(&temp, noisePointer, &(*pair));

	  point_t fi;
	  point_t fj;

	  // compute the particle-expressions
	  m_1stparticleFactor(&fi, noisePointer, &(*pair));
	  m_2ndparticleFactor(&fj, noisePointer, &(*pair));

	  // loop necessary because operator* of math_vector_t does scalar product
	  for(size_t i = 0; i < SPACE_DIMS; ++i) {
	    fi[i] *= temp[i];
	    fj[i] *= temp[i];
	  }

	  double sqrtDt = sqrt(m_controller->dt());
	  fi /= sqrtDt;
	  fj /= sqrtDt;

	  Particle* first = pair->firstPart(); 
	  Particle* second = pair->secondPart();
	  //             MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "temp = " << temp);            
	  

/* 	  if (first->mySlot==100) MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "0first: tag-data BEFORE = " << first->tag.pointByOffset(m_slots.first) << ", m_symmetry = " << m_symmetry); */
/* 	  if (second->mySlot==100) MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "0second: tag-data BEFORE = " << second->tag.pointByOffset(m_slots.second) << ", m_symmetry = " << m_symmetry); */

          
	  if(pair->actsOnFirst()) {
	      
#ifndef _OPENMP
          
/*  	      MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "tag-data BEFORE = " << first->tag.pointByOffset(m_slots.first)); */


	      first->tag.pointByOffset(m_slots.first) += fi;

/* 	      MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "tag-data AFTER = " << first->tag.pointByOffset(m_slots.first));             */


#else
	      // FIXME: parallelise!
              first->tag.pointByOffset(m_slots.first) += fi;
#endif
	      
	      //     MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "AFTER: first->tensor = " << first->tag.pointByOffset(m_slots.first));            

	  }
	  
	  if(pair->actsOnSecond()) {
	    
#ifndef _OPENMP
	    second->tag.pointByOffset(m_slots.second) += m_symmetry*fj;
#else
	    // FIXME: parallelise
	    second->tag.pointByOffset(m_slots.second) += m_symmetry*fj;
#endif
	      
	    //     MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "AFTER: first->tensor = " << first->tag.pointByOffset(m_slots.first));

	  }


/* 	  if (first->mySlot==100) MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "0first: tag-data AFTER = " << first->tag.pointByOffset(m_slots.first)); */
/* 	  if (second->mySlot==100) MSG_DEBUG("BondedPairParticleTensorNoiseVector::compute", "0second: tag-data AFTER = " << second->tag.pointByOffset(m_slots.second)); */
	  
        }

	/*!
         * Returns the symbol name as defined in the input file.
	 */
        virtual string myName() {
          return m_symbolName;
        }

	/*!
         * Returns the name of the bonded list this \a ValCalculator computes for.
	 */
        virtual string listName() {
          return m_listName;
        }

	/*!
         * Register the computed symbol
	 */
        virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) 
        {
          throw gError("BondedPairParticleTensorNoiseVector::setSlots", "should not have been called! Contact the programmer.");
        }

	/*!
         * Setup this Calculator
	 */
        virtual void setup();    

	/*!
	 * copies the members of this class to \a vc
	 */
	virtual void copyMembersTo(ValCalculator* vc);



};

#endif
