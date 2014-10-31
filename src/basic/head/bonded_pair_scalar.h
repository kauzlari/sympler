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


#ifndef __BONDED_PAIR_SCALAR_H_
#define __BONDED_PAIR_SCALAR_H_

#include "bonded_pair_arbitrary.h"
//#include "function_pair.h"
#include "general.h"

using namespace std;

/*! Module for calculation of a scalar Symbol stored in a pair*/

class BondedPairScalar : public BondedPairArbitrary {
  
 protected:
  
  virtual ValCalculatorPair* copyMySelf() {
    return new BondedPairScalar(*this);
  }

  /*!
   * Initialise the property list
   */
  virtual void init();

 public:
  
  /*!
   * Constructor for Node hierarchy//
   */
  BondedPairScalar(/*Node*/Simulation* parent);
  
  /*!
   * Destructor//
   */
  virtual ~BondedPairScalar();

  /*!
   * Setup this Calculator
   */
  virtual void setup();
  
  /*!
   * Returns the symbol name as defined in the input file.
   */
  virtual string myName() {
    return m_symbolName;
  }

// Now (2014-10-31) in BondedPairArbitrary
//  /*! 
//   * set the slots for the pairs
//   */
//  virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp);

  void mySlot(size_t& slot) const {
    slot = m_slot;
  }
  
  /*! 
   * Compute the user defined expression for pair dis
   */
#ifdef _OPENMP
  virtual void compute(Pairdist* dis, int threadNum) {
    // does the same as in serial mode at the moment; 
    // it writes into pairs, so should work
    compute(dis);
  }
#endif
  virtual void compute(Pairdist* dis) {
      
/*       double temp; */
      
      // compute the expression for the pair function
/*       m_function(&temp, dis); */
      m_function(&(dis->tag.doubleByOffset(m_slot)), dis);
/*       dis->tag.doubleByOffset(m_slot)=temp; */

      // MSG_DEBUG("BondedPairScalar::compute", "In pair with value " << temp);
    
  }

};

#endif 
