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



#ifndef __VAL_CALCULATOR_PART_H
#define __VAL_CALCULATOR_PART_H

#include "symbol.h"

// #define VC_MAX_STAGE 1

class Pairdist;
class ColourPair;


/*!
 * Base class for all \a ValCalculator s that store information within a \a Particle
 */
class ValCalculatorPart : public ValCalculator
{
  public:

  /*!
   * Constructor
   */
    ValCalculatorPart(string symbol)
  : ValCalculator(symbol)
    {
//       m_symbolName = symbol;
    }

  /*!
     * Constructor for Node hierarchy
   */
    ValCalculatorPart(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorPart() {
    }

  /*!
     * Will throw an exception
   */
    virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp)
    {
      throw gError("ValCalculatorPart::setSlot(ColourPair*, size_t&) should "
          "not be called!");
//  /*!
//     * Setup the Calculator if used in a Node hierarchy
      //   */
//    virtual void setup();
    }

  /*!
     * Will throw an exception
   */
    void mySlot(size_t& slot) const
    {
      throw gError("ValCalculatorPart::mySlot(size_t&) should not be called!");
    }

  /*!
   * Return the offset for the data stored by this calculator in the \a Particle s
   */
    void mySlots(pair<size_t, size_t> &theSlots)
    {
      theSlots = m_slots;
    }


  /*!
   * Return the offset for the data stored by this calculator in the \a Particle s
   */
    virtual pair<size_t, size_t> vcSlots()
    {
      return m_slots;
    }

#ifdef _OPENMP
  /*!
   * Return the offset for the data stored by this calculator in the \a Particle s
   */
    virtual vector<pair<int, int> > &copySlots() {
      return m_copy_slots;
    }

  /*!
   * Return the offset for the data inside the copy vector in the \a Particle s tag.
   */
    virtual pair<int, int> &vectorSlots() {
      return m_vector_slots;
    }

  /*!
   * Merging the copies among different threads (processors) together
   */
    virtual void mergeCopies(ColourPair* cp, int thread_no) = 0;
#endif

   /*!
    * Setup this Calculator
    */
    virtual void setup();

  protected:
  /*!
   * The tag offset of the data stored in a \a Particle
   */
    pair<size_t, size_t> m_slots;

#ifdef _OPENMP
  /*!
   * The tag offset of the copy data stored in a \a Particle.
   * The vector size equals the number of threads used.
   */
    vector<pair<int, int> > m_copy_slots;

    /*!
     * The offset to the copy-slots inside a vector
     */
    pair<int, int> m_vector_slots;
#endif

  /*!
   * Initialise the property list
   */
    virtual void init();
};


#endif
