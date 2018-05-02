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


#ifndef __VAL_CALCULATOR_NEG_DKERNEL_DIVR_H
#define __VAL_CALCULATOR_NEG_DKERNEL_DIVR_H

#include "val_calculator.h"
#include "weighting_function.h"

//---- newclass -----------------------------------------------------------

/*!
 * Saves -(1/r)*(dW/dr) for each pair of particles, where W is an interpolation function (the kernel) and r is the distance
 */
class ValCalculatorNegDKernelDivr : public ValCalculatorPair
{
  protected:
    /*!
     * The weighting function to be used for the local density calculation
    */
    WeightingFunction *m_wf;

    /*!
    * Name of the weighting function
    */
    string m_wfName;

    /*!
    * Initialise the property list
    */
    virtual void init();

    /*!
     * Helper function for returning a copy of itself
    */
    virtual ValCalculatorPair* copyMySelf()
    {
      return new ValCalculatorNegDKernelDivr(*this);
    }

  public:
  /*!
   * Constructor
   * @param wf Weighting function used for computation
   * @param symbol The symbol name
  */
    ValCalculatorNegDKernelDivr(WeightingFunction *wf, string symbol);

  /*!
     * Constructor for Node hierarchy
  */
    ValCalculatorNegDKernelDivr(/*Node*/Simulation* parent);

  /*!
     * Destructor
  */
    virtual ~ValCalculatorNegDKernelDivr() {
}

    /*!
     * Register this calculator and save the offset for the data stored
     * in a \a Pairdist in the argument \a slot
    */
    void setSlot(ColourPair* cp, size_t& slot, bool oneProp);

    /*!
     * Compute the kernel function for the \a Pairdist \a pD
     * @param pD \a Pairdist for which to compute the kernel function
     */
#ifndef _OPENMP
    virtual void compute(Pairdist* pD)
#else
    virtual void compute(Pairdist* pD, int thread_no)
#endif
    {
      point_t dummy = {{{0, 0, 0}}};
      pD->tag.doubleByOffset(m_slot) = m_wf->weight(pD, dummy, NULL);
    }

  /*!
     * Returns "kernel"
  */
    virtual string myName() {
      return "negDkernelDivr";
}

    /*!
     * Setup this Calculator
    */
    virtual void setup();

};

#endif
