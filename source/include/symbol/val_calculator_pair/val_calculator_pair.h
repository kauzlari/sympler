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



#ifndef __VAL_CALCULATOR_PAIR_H
#define __VAL_CALCULATOR_PAIR_H

#include "symbol.h"

// #define VC_MAX_STAGE 1

class Pairdist;
class ColourPair;


/*!
 * Base class for all \a ValCcalculator s that store information within a \a Pairdist
 */
class ValCalculatorPair : public ValCalculator
{

 protected:

/*!
 * Should all colour combinations be considered. This disables \a m_species.
 */
  bool m_allPairs;


 public:

    /*!
 * Constructor
     */
  ValCalculatorPair(string symbol)
  : ValCalculator(symbol)
  {}

    /*!
   * Constructor for Node hierarchy
     */
  ValCalculatorPair(/*Node*/Simulation* parent);

  /*!
   * Destructor
   */
  virtual ~ValCalculatorPair() {
  }

  virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp);

  virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) {
    throw gError
      ("ValCalculatorPair::setSlots", "(ColourPair*, "
       "pair<size_t, size_t>&, bool ) should not be called! Contact the programmers.");
  }

  void mySlot(size_t& slot) const
  {
    slot = m_slot;
  }

  void mySlots(pair<size_t, size_t> &theSlots)
  {
    throw gError("ValCalculatorPair::mySlots(pair<size_t, size_t>*) should "
                 "not be called!");
  }

  /*!
  * Setup this Calculator
  */
  virtual void setup();

//  /*!
//  * Setup the Calculator if used in a Node hierarchy
  //  */
//  virtual void setup();

protected:

  /*!
   * The tag offset of the data stored in a \a Pairdist
   */
  size_t m_slot;

  /*!
   * Initialise the property list
   */
  virtual void init();

  /*!
   * Helper function for returning a copy of itself
   */
  virtual ValCalculatorPair* copyMySelf() = 0;
};


#endif
