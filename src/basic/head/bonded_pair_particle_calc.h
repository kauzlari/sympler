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


#ifndef __BONDED_PAIR_PARTICLE_CALC_H
#define __BONDED_PAIR_PARTICLE_CALC_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"


/*!
 * Parent class for computing properties for the 
 * particles, which need summation over bonded pairs. 
 */
class BondedPairParticleCalc : public ValCalculatorPart
{
 protected:
  
  /*!
   * name of the bonded list, this calculator belongs to
   */
  string m_listName;
  
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
   * Helper function for polymorphic copying
   */
  virtual ValCalculator* copyMySelf() = 0;
  
  /*!
   * copies the members of this class to \a vc
   */
  virtual void copyMembersTo(ValCalculator* vc);

  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * There are subclasses, which indeed do not leave the list in its original state
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
  {

  }

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
  BondedPairParticleCalc(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~BondedPairParticleCalc();

#ifdef _OPENMP
  /*!
   * Merge the copies of all threads together
   */
  virtual void mergeCopies(ColourPair* cp, int thread_no);
#endif
  
  /*!
   * Compute the user defined expression for pair \a pD
   * @param pD \a Pairdist whose contribution we calculate
   */
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
   * Register the computed symbol
   */
  virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) 
  {
    throw gError("BondedPairParticleVector::setSlots", "should not have been called! Contact the programmer.");
  }
  
  /*!
   * Setup this Calculator
   */
  virtual void setup();
    
};

#endif
