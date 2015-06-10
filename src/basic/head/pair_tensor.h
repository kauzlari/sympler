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



#ifndef __PAIR_TENSOR_H_
#define __PAIR_TENSOR_H_


#include "pair_arbitrary.h"

using namespace std;

/*! Module for calculation of a square matrix Symbol which is stored in each pair (within the cutoff).*/

class PairTensor: public PairArbitrary {

 protected:  
  
  virtual ValCalculatorPair* copyMySelf() {
    return new PairTensor(*this);
  }
  /*!
   * Initialise the property list
   */
  virtual void init();
  
 public:
  
  PairTensor(/*WeightingFunction *wf,*/string symbol);
  /*!
   * Constructor for Node hierarchy//
   */
  PairTensor(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~PairTensor();
  
  /*!
   * Setup this Calculator
   */
  virtual void setup();
  
  /*! Compute the user defined expression for pair \a dis.*/
  virtual void compute(Pairdist* dis) {
    
    if (dis->abs() < m_cutoff) {
  
      tensor_t& tensor = dis->tag.tensorByOffset(m_slot);
          
      // compute the expression for the pair function
      m_function(&tensor, dis);
    }
  }

#ifdef _OPENMP
  /*! Compute the user defined expression for pair \a dis -- parallel version. Should be that simple because only looping over and writing into pairs (not particles of the pair!)*/
  virtual void compute(Pairdist* dis, int threadNum) {
    // does the same as in serial mode at the moment; 
    // it writes into pairs, so should work
    compute(dis);
  }
#endif


};

#endif
