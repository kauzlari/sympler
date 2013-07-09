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


#ifndef __VAL_CALCULATOR_DIRICHLETBC_H
#define __VAL_CALCULATOR_DIRICHLETBC_H

#include "val_calculator_BC.h"

// #include "wall.h"

/* class Wall; */

//---- newclass -----------------------------------------------------------

/*!
 * Saves the pair-specific value of a quantity of the boundary particle used
 * for applying a Dirichlet boundary condition (BC) in each pair of particles.
 * Which kind of quantity we are dealing with is specified by the child
 * classes.
 * Notice: Right now (2007/10/24) This class only implements the velocity
 * Dirichlet-BC
 */
class ValCalculatorDirichletBC : public ValCalculatorBC
{
  protected:

    /*!
     * Initialise the property list
     */
    virtual void init();

    /*!
     * Helper function for returning a copy of itself
     */
    virtual ValCalculatorPair* copyMySelf()
    {
      return new ValCalculatorDirichletBC(*this);
    }

  public:
  /*!
   * Constructor
   * @param symbol The symbol name
   */
    ValCalculatorDirichletBC(string symbol);

  /*!
     * Constructor for Node hierarchy
   */
    ValCalculatorDirichletBC(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorDirichletBC() {
    }

    /*!
     * Register this calculator and save the offset for the data stored
     * in a \a Pairdist in the argument \a slot
     */
/*     void setSlot(ColourPair* cp, size_t& slot, bool oneProp); */

    /*!
     * Compute the boundary value for the \a Pairdist \a pD
     * @param pD \a Pairdist for which to compute the boundary value
     */
#ifndef _OPENMP
    virtual void compute(Pairdist* pair); //FIXME: inline ?
#else
    virtual void compute(Pairdist* pD, int thread_no);
#endif

  /*!
     * Returns "dirichletBC"
   */
    virtual string myName() {
      return "dirichletBC";
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();

    /*!
     * For each wall-particle, find the walls in range
     */
/*     virtual void setupAfterParticleCreation(); */
};

#endif
