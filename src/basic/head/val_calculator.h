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



#ifndef __VAL_CALCULATOR_H
#define __VAL_CALCULATOR_H

#include "symbol.h"

/* This file contains definitions of

ValCalculator
ValCalculatorPair
ValCalculatorPart
ValCalcPartArbitrary

*/

class Pairdist;
class ColourPair;

/*!
 * Calculator for cached properties computed during a loop over pairs
 */
class ValCalculator : public Symbol
{
  
 protected:
  
  /*!
   * The species of the ColourPair this Calculator should belong to
   */
  pair<string, string> m_species;

#ifdef _OPENMP
  /*!
   * Is this a particle or a pair calculator.
   */
  bool m_particleCalculator;
#endif
  
  /*!
   * Initialise the property list
   */
  virtual void init();

  /*!
   * Helper function which removes indices and brackets from single 
   * terms in \a Function expressions. E.g.: "[rij]" becomes "r" or 
   * "rhoi" becomes "rho"
   * @param[out] name Term from a \a Function expression to be cleaned
   */
  virtual void cleanSymbol(string& name) const;

  
 public:
  
  /*!
   * Constructor
   */
  ValCalculator(string symbol);

  /*!
   * Constructor
   */
  ValCalculator(/*Node*/Simulation* parent);

  /*!
   * Destructor
   */
  virtual ~ValCalculator() {
  }

  /*!
   * Setup the variables of this \a Symbol
   */
  virtual void setup();
  
  /*!
   * Compute cached properties for \a Pairdist \a pD
   * @param pD The \a Pairdist for which to compute the cached properties
   */
#ifndef _OPENMP
  virtual void compute(Pairdist* pD) = 0;
#else
  virtual void compute(Pairdist* pD, int thread_no) = 0;
  
  virtual bool particleCalculator() {
    return m_particleCalculator;
  }
#endif
  
  /*!
   * Return a string identifier for this calculator
   */
  virtual string myName() = 0;

  /*!
   * Return the offset for the data stored by this calculator in a \a Pairdist
   */
  virtual void mySlot(size_t& slot) const = 0;
  
  /*!
   * Return the offset for the data stored by this calculator in the 
   * \a Particle s
   */
  virtual void mySlots(pair<size_t, size_t> &theSlots) = 0;

  /*!
   * Register this calculator and save the offset for the data stored
   * in a \a Pairdist in the argument \a slot
   */
  virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp) = 0;

  /*!
   * Register this calculator and return the offset for the data stored
   * in the \a Particle s
   */
  virtual void setSlots
    (ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) = 0;
    
  /*!
   * Return cutoff.
   * @return Cutoff returned here serves for child-debugging only
   */
  virtual double cutoff() {
    return -1;
  }

  /*!
   * Return the first species this calculator writes into. Here, the 
   * default behaviour is implemented
   */
  virtual string firstWriteSpecies() {
    return m_species.first;
  }

  /*!
   * Return the second species this calculator writes into. Here, the default behaviour is implemented
   */
  virtual string secondWriteSpecies() {
    return m_species.second;
  }

};


/*!
 * Base class for all \a ValCcalculator s that store information within 
 * a \a Pairdist
 */
class ValCalculatorPair : public ValCalculator
{

 protected:

  /*!
   * Should all colour combinations be considered. This disables \a m_species.
   */
  bool m_allPairs;

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

  /*!
   * Helper function for checking that a copy of \a this 
   * \a ValCalculatorPair is a safe one. This check includes here by 
   * default only assertions. But subclasses could do more, for example 
   * also register the copy to some computations (precompute, 
   * setupBeforeParticleCreation, or similar)
   * @param[in] vc Copy of a \a ValCalculatorPair to check
   */
  virtual void makeCopySafe(ValCalculatorPair* vc) const;

  
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
  
  virtual void setSlots
    (ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) {
    throw gError
      ("ValCalculatorPair::setSlots for module " + className(),
       "(ColourPair*, pair<size_t, size_t>&, bool ) should not be "
       "called! Contact the programmer of the module.");
  }

  void mySlot(size_t& slot) const
  {
    slot = m_slot;
  }
  
  void mySlots(pair<size_t, size_t> &theSlots)
  {
    throw gError("ValCalculatorPair::mySlots for module " + className(),
		 "should not be called!");
  }

  /*!
   * Setup this Calculator if used in a Node hierarchy
   */
  virtual void setup();

};


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
    {}
  
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
  virtual void setSlot(ColourPair* cp, size_t& slot, bool oneProp) {
    throw gError("ValCalculatorPart::setSlot for module " + className(),
		 "should not be called!");
  }

  /*!
   * Will throw an exception
   */
  void mySlot(size_t& slot) const
  {
    throw gError("ValCalculatorPart::mySlot for module " + className(),
		 "should not be called!");
  }

  /*!
   * Return the offset for the data stored by this calculator in the 
   * \a Particle s
   * @param[out] theSlots Storage for the returned offsets (may be 
   * different for different colours)
   */
  void mySlots(pair<size_t, size_t> &theSlots)
  {
    theSlots = m_slots;
  }
  
  /*!
   * Return the offset for the data stored by this calculator in the 
   \a Particle s
   * @return The offsets for both particles (may be different for 
   * different colours)
   */
  virtual pair<size_t, size_t> vcSlots() {
    return m_slots;
  }

#ifdef _OPENMP
  /*!
   * Return the offset in the \a Particle tag to the data copies for 
   * each thread stored by this calculator 
   */
  virtual vector<pair<int, int> > &copySlots() {
    return m_copy_slots;
  }

  /*!
   * Return the actual slot to the data inside the copy vector 
   * (accessed in the \a Particle tag by \a m_copy_slots)
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
     * The offset in the \a Particle tag to the data copies for 
     * each thread stored by this calculator 
     */
    vector<pair<int, int> > m_copy_slots;
    
    /*!
     * The actual slot to the data inside the copy vector 
     * (accessed in the \a Particle tag by \a m_copy_slots)
     */
    pair<int, int> m_vector_slots;
    
#endif
    
    /*!
     * Initialise the property list
     */
    virtual void init();
    
};

/*!
 * Parent class of all \a ValcalculatorPart which loop over the 
 * non-bonded neighbour list
 */
class NonBondedPairParticleCalculator : public ValCalculatorPart
{
  
 protected:

  /*!
   * Should all colour combinations be considered. This disables \a m_species.
   */
  bool m_allPairs;

  /*!
   * initialise this Caluclator
   */
  virtual void init();

  
 public:
  
  /*!
   * Constructor
   */
  NonBondedPairParticleCalculator(string symbol);
  
  
  /*!
   * Constructor for Node hierarchy
   */
  NonBondedPairParticleCalculator(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~NonBondedPairParticleCalculator() {
  }

};


#endif
