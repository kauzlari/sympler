/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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


#ifndef __PAIR_RAND_VECTOR_H_
#define __PAIR_RAND_VECTOR_H_

#include "pair_rand_arbitrary.h"

using namespace std;

/*! Module for calculation of a random vector Symbol [pairFactor]*[ranVector] stored per non-bonded pair of particles, where [pairFactor] is a user-defined runtime-compiled vector expression and [ranVector] is a vector of independent Gaussian random numbers with zero mean and unit variance. Note that you can modify the mean and variance through a suitable [pairFactor]. [pairFactor] and [ranVector] are multiplied componentwise.*/

class PairRandVector : public PairRandArbitrary {

 protected:
  
  /*!
   * For copying an instance of this class
   */  
  virtual ValCalculatorPair* copyMySelf() {
    return new PairRandVector(*this);
  }
  /*!
   * Initialise the property list
   */
  virtual void init();

 public:
  
  /*!
   * Old, usually not used constructor
   */
  PairRandVector(string symbol);

  /*!
   * Constructor for Node hierarchy
   */
  PairRandVector(/*Node*/Simulation* parent);
  
  /*!
   * Destructor//
   */
  virtual ~PairRandVector();
    
  
  /*Compute the user defined expression for pair pD*/
#ifdef _OPENMP
  virtual void compute(Pairdist* dis, int threadNum) {
    // does the same as in serial mode at the moment; 
    // it writes into pairs, so should work
    compute(dis);
  }
#endif
  virtual void compute(Pairdist* dis) {
    if (dis->abs() < m_cutoff) {

      point_t& vector =  dis->tag.pointByOffset(m_slot);
      
      // compute the expression for the pair function
      m_function(&vector, dis);

      // multiply Gaussian random numbers
      for(size_t i = 0; i < SPACE_DIMS; ++i)
	vector[i] *=  m_rng.normal(1);

    }    
  }
  
};

#endif 
