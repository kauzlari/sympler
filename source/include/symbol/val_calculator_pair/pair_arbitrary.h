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



#ifndef __PAIR_ARBITRARY_H_
#define __PAIR_ARBITRARY_H_


#include "val_calculator.h"
#include "val_calculator_pair.h"
#include "function_pair.h"
#include "general.h"

using namespace std;

/*! Module for calculation of a square matrix Symbol which is stored in each pair (within the cutoff).*/

class PairArbitrary: public ValCalculatorPair {

 protected:

  /*!
   * Cut-off radius for the pair summation
   */
  double m_cutoff;

  /*!
   * A string for m_function;
   */
  string m_expression_string;
  
  /*!
   * The function pair computing the user defined pair factor
   */
  FunctionPair m_function;
  
  /*!
   * This string holds the symbols, which are not waited for to be computed beforehand
   */
  /* string m_oldSymbols; */
  
  
  virtual ValCalculatorPair* copyMySelf() = 0;

  /*!
   * Initialise the property list
   */
  virtual void init();

  /*!
   * The returned string contains those terms from runtime compiled expressions, 
   * which should be ignored when determining the stage. The expressions are separated by " | ".
   * An "empty" string must have the form "---".
   */
  /* virtual string usedSymbolsIgnoredForStaging() const { */
  /*   return m_oldSymbols; */
  /* } */

  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols);

  /*!
   * Returns the strings of those \a Symbols that the given class depends on
   * due to hard-coded reasons (not due to runtime compiled expressions).
   * @param usedSymbols List to add the strings to.
   */
  virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
  {
    
  }

  
 public:
  
  /*!
   * Old-style constructor. Most likely not used anymore and to be removed in the future.
   */
  PairArbitrary(string symbol);

  /*!
   * Constructor for Node hierarchy
   */
  PairArbitrary(/*Node*/Simulation* parent);
  
  /*!
   * Destructor//
   */
  virtual ~PairArbitrary();
  
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
  
  /*! set the slots for the pairs */
  virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp);
  
  /*! Slots for each particle; we need a slot per pair, so this should not be called */
  virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots,
			bool oneProp) {
    throw gError
      ("PairArbitrary::setSlots(ColourPair*, "
       "pair<size_t, size_t>&, bool) for " + className(), "should not be called! Contact the programmers.");
  }
  
  /*! set the slots for the pairs */
  void mySlot(size_t& slot) const {
    slot = m_slot;
  }
  
  /*! Slots for each particle; we need a slot per pair, so this should not be called */
  virtual void mySlots(pair<size_t, size_t> &theSlots) {
    throw gError("PairArbitrary::mySlots(pair<size_t, size_t>*) for " + className(),"should "
		 "not be called! Contact the programmers.");
  }
  
};


#endif
