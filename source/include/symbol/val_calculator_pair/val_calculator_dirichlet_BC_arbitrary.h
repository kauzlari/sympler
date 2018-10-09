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

#ifndef __VAL_CALCULATOR_DIRICHLET_BC_ARBITRARY_H
#define __VAL_CALCULATOR_DIRICHLET_BC_ARBITRARY_H

#include "val_calculator_BC.h"

#include "particle.h"
#include "pairdist.h"

/*!
 * Saves the pair-specific value of an arbitrary variable of the boundary
 * particle used for applying a Dirichlet boundary condition (BC) in each pair
 * of particles.
 */
class ValCalculatorDirichletBCArbitrary : public ValCalculatorBC
{

	protected:

		/*!
		 * Name of the scalar the BC is computed for
		 */
		string m_varName;

		/*!
		 * Memory offset to the scalar descibed by \a m_scalarName
		 */
		pair<size_t, size_t> m_varOffset;
   
    /*!
     * Initialise the property list
     */
    virtual void init();
  
    /*!
     * Helper function for returning a copy of itself
     */
    virtual ValCalculatorPair* copyMySelf() = 0;

    /*!
     * Returns the strings of those \a Symbols that the given class depends on
     * due to hard-coded reasons (not due to runtime compiled expressions).
     * @param usedSymbols List to add the strings to.
     */
    virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
    {
      usedSymbols.push_back(m_varName);
    }

  
  public:
  /*!
   * Constructor
   * @param symbol The symbol name
   */
    ValCalculatorDirichletBCArbitrary(string symbol);

  /*!
     * Constructor for Node hierarchy
   */
    ValCalculatorDirichletBCArbitrary(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorDirichletBCArbitrary() {
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
    virtual void compute(Pairdist* pD) = 0;
#else
      virtual void compute(Pairdist* pD, int thread_no) = 0;
#endif


  /*!
     * Returns an identifier
   */
    virtual string myName() = 0;

    /*!
     * Setup this Calculator
     */
    virtual void setup();

};

#endif
