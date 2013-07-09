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



#ifndef __VAL_CALCULATOR_R_I6_H
#define __VAL_CALCULATOR_R_I6_H

#include "val_calculator.h"

/*!
 * Return the reciprocal distance to the power of six (for the Lennard-Jones potential)
 */
class ValCalculatorRi6 : public ValCalculatorPair
{
  protected:

  /*!
   * Initialise the property list
   */
    virtual void init()
    {
      m_properties.setClassName("ValCalculatorRi6");

      m_properties.setDescription("Saves 1/r^6, i.e., the 6th power of the inverse of the distance in each pair of particles.");

      STRINGPC
          (symbol, m_symbolName,
           "Name of the symbol for the calculated value.");

      #ifdef _OPENMP
      m_particleCalculator = false;
      #endif

    }

    /*!
     * Helper function for returning a copy of itself
     */
    virtual ValCalculatorPair* copyMySelf()
    {
      return new ValCalculatorRi6(*this);
    }

  public:
  /*!
   * Constructor
   */
  ValCalculatorRi6()
  : ValCalculatorPair("ri6")
  {
    m_stage = 0;
  }

  /*!
   * Constructor for \a Node hierarchy
     */
  ValCalculatorRi6(/*Node*/Simulation* parent)
  : ValCalculatorPair(parent)
  {
    m_stage = 0;
    init();
  }

  /*!
   * Destructor
   */
  virtual ~ValCalculatorRi6() {
  }

  /*!
   * Compute 1/r^6
   */
#ifndef _OPENMP
  virtual void compute(Pairdist* pD)
#else
  virtual void compute(Pairdist* pD, int thread_no)
#endif
  {
    double temp = 1/pD->absSquare();

    pD->tag.doubleByOffset(m_slot) = temp*temp*temp;
  }

  /*!
   * Return "ri6"
   */
  virtual string myName() {
    return "ri6";
  }

  /*!
   * Setup this Calculator
     */
  virtual void setup()
  {
    ValCalculatorPair::setup();
  }

};

#endif
