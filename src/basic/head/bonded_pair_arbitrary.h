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


#ifndef __BONDED_PAIR_ARBITRARY_H
#define __BONDED_PAIR_ARBITRARY_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "function_pair.h"


/*!
 * Parent class for computing properties for the bonded pairs. 
 */
class BondedPairArbitrary : public ValCalculatorPair
{
  protected:

  /*!
   * name of the bonded list, this calculator belongs to
   */
  string m_listName;

  /*!
   * The mathematical pair-expression to be computed
   */
  string m_expression;
  
  /*!
   * The \a FunctionPair computing the user defined pair factor
   */
  FunctionPair m_function;
  
  /*!
   * This string holds the symbols, which are not waited for to be computed beforehand
   */
  string m_oldSymbols;

  /*!
   * Initialise the property list
   */
  virtual void init();
  

  /*!
   * Helper function for polymorphic copying
   */
  virtual ValCalculatorPair* copyMySelf() = 0;
  /*         { */
  /*           return new BondedPairArbitrary(*this); */
  /*         } */
  
    /*!
     * copies the members of this class to \a vc
     */
    virtual void copyMembersTo(ValCalculator* vc);

  public:
    
    /*!
     * Constructor for the \a Node hierarchy
     */
    BondedPairArbitrary(/*Node*/Simulation* parent);
    
    /*!
     * Destructor
     */
    virtual ~BondedPairArbitrary();
  
      /*!
       * Compute the user defined expression for pair \a pD
       * @param pD \a Pairdist whose contribution we calculate
       */

	//set the slots for the pairs
	virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp);

	//Slots for each particle,we need a slot per pair,so this should not be called
	virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots,
			      bool oneProp) {
	  throw gError
	    ("BondedPairArbitrary::setSlots", "(ColourPair*, "
	     "pair<size_t, size_t>&, bool ) should not be called! Contact the programmers.");
	}

	void mySlot(size_t& slot) const {
		slot = m_slot;
	}

	virtual void mySlots(pair<size_t, size_t> &theSlots) {
		throw gError("BondedPairArbitrary::mySlots(pair<size_t, size_t>*) should "
				"not be called! This is a fatal error to be reported to a programmer.");
	}

#ifndef _OPENMP

	virtual void compute(Pairdist* pD) = 0;

#else

	virtual void compute(Pairdist* pD, int thread_no) = 0;

#endif
    
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
     * Setup this Calculator
     */
    virtual void setup();

    /*!
     * Diffenrently to the function in \a Symbol, this class really has
     * to determine its stage during run-time
     */
    virtual bool findStage();
    
    /*!
     * Diffenrently to the function in \a Symbol, this class really has
     * to determine its stage during run-time
     */
    virtual bool findStage_0();
    
};

#endif
