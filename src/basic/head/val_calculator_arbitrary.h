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


#ifndef __VAL_CALCULATOR_ARBITRARY_H
#define __VAL_CALCULATOR_ARBITRARY_H

// #include "general.h"
// #include "simulation.h"
// #include "manager_cell.h"
#include "val_calculator.h"
#include "function_pair.h"
#include "colour_pair.h"


/*!
 * Parent class with functions to compute completely user-defined properties for the
 * particles, which need pair summation
 */
class ValCalculatorArbitrary : public NonBondedPairParticleCalculator
{
  protected:

  /*!
   * Cut-off radius for the pair summation
   */
    double m_cutoff;

   /*!
    * The mathematical pair-expression to be computed
        */
    string m_expression;

    /*!
     * The \a FunctionPair computing the user defined pair factor
     */
    FunctionPair m_function;

    /*!
     * An additional factor for the 1st particle
     */
    FunctionPair m_1stparticleFactor;

    /*!
     * An additional factor for the 2nd particle
     */
    FunctionPair m_2ndparticleFactor;

   /*!
     * The mathematical expression for \a m_1stparticleFactor
    */
    string m_1stPExpression;

   /*!
     * The mathematical expression for \a m_2ndparticleFactor
    */
    string m_2ndPExpression;

    /*!
    * This string holds the symbols, which are not waited for to be computed beforehand
    */
    string m_oldSymbols;

    /*! This is now a member of \a Symbol
     * Is this calculator allowed to overwrite already existing symbols
     * with name \a m_symbolName ?
     */
/*     bool m_overwrite; */

    /*!
    * How does the pair expression behave under index interchange?
    * 1 means symmetric behaviour
    * -1 means antisymmetric behaviour
    */
    int m_symmetry;

  /*!
     * Initialise the property list
   */
    virtual void init();

     /*!
     * Return a copy of the current object
      */
    virtual ValCalculator* copyMySelf() = 0;

    /*!
     * The returned string contains those terms from runtime compiled expressions, 
     * which should be ignored when determining the stage. The expressions are separated by " | ".
     * An "empty" string must have the form "---".
     */
    virtual string usedSymbolsIgnoredForStaging() const {
      return m_oldSymbols;
    }

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
     * Constructor for the \a Node hierarchy
     */
    ValCalculatorArbitrary(/*Node*/Simulation* parent);
    
    /*!
     * Destructor
     */
    virtual ~ValCalculatorArbitrary() {
    }
    
    /*!
     * Setup this Calculator
     */
    virtual void setup();
    
    /*!
     * Compute the user defined expression for pair \a pD
     * @param pD \a Pairdist whose contribution is calculated
     */
#ifndef _OPENMP
    virtual void compute(Pairdist* pD) = 0;
#else
    virtual void compute(Pairdist* pD, int thread_no) = 0;
#endif
    
    
#ifdef _OPENMP
    /*!
     * Merge the copies of all threads together
     */
    virtual void mergeCopies(ColourPair* cp, int thread_no) = 0;
#endif
    
    /*!
     * Returns the symbol name as defined in the input file.
     */
    virtual string myName() {
      return m_symbolName;
    }
    
    /*!
     * Register all degrees of freedom
     */
    virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp)
    {
      throw gError("ValcalculatorArbitrary::setSlots", "should not have been called! Contact the programmer.");
    }
    
};

#endif
