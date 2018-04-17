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



#ifndef __VAL_CALCULATOR_H
#define __VAL_CALCULATOR_H

#include "symbol.h"

// #define VC_MAX_STAGE 1

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

  // now moved to NonBondedPairParticleCalculator and ValCalculatorPair

  //  /*!
  //  * Should all colour combinations be considered. This disables \a m_species.
  //  */
  //  bool m_allPairs;

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
   * Helper function which removes indices and brackets from single terms in \a Function 
   * expressions. E.g.: "[rij]" becomes "r" or "rhoi" becomes "rho"
   * @param name Single term from a \a Function expression
   */
  virtual void cleanSymbol(string& name) const;

  
public:
  /*!
 * Constructor
   */
  ValCalculator(string symbol);/*: m_stage(0) {
}*/

  /*!
   * Constructor
   */
  ValCalculator(/*Node*/Simulation* parent);/*: m_stage(0) {
}*/

  /*!
   * Destructor
   */
  virtual ~ValCalculator() {
  }

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
   * Return the offset for the data stored by this calculator in the \a Particle s
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
  virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool
			oneProp) = 0;
  
  //  /*!
  //   * Return the stage at which this calculator is being called.
  //   */
  //  virtual size_t stage() {
  //    return m_stage;
  //  }
  
  /*!
   * return cutoff
   */
  virtual double cutoff() {
    return -1;
  }

  /*!
   * Return the first species this calculator writes into. Here, the default behaviour is implemented
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
 * Base class for all \a ValCcalculator s that store information within a \a Pairdist
 */
class ValCalculatorPair : public ValCalculator
{

 protected:

/*!
 * Should all colour combinations be considered. This disables \a m_species.
 */
  bool m_allPairs;


 public:

    /*!
 * Constructor
     */
  ValCalculatorPair(string symbol)
  : ValCalculator(symbol)/*, m_persistency(false)*/
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

  virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) {
    throw gError
      ("ValCalculatorPair::setSlots", "(ColourPair*, "
       "pair<size_t, size_t>&, bool ) should not be called! Contact the programmers.");
  }

  void mySlot(size_t& slot) const
  {
    slot = m_slot;
  }

  void mySlots(pair<size_t, size_t> &theSlots)
  {
    throw gError("ValCalculatorPair::mySlots(pair<size_t, size_t>*) should "
                 "not be called!");
  }

  /*!
  * Setup this Calculator
  */
  virtual void setup();

//  /*!
//  * Setup the Calculator if used in a Node hierarchy
  //  */
//  virtual void setup();

protected:
  /*!
   * The tag offset of the data stored in a \a Pairdist
   */
  size_t m_slot;

  // see CONVENTION5 for rule about persistencies

  /*
  * If this is true, the symbol cannot be cleared
  */
//   bool m_persistency;

  /*!
   * Initialise the property list
   */
  virtual void init();

  /*!
   * Helper function for returning a copy of itself
   */
  virtual ValCalculatorPair* copyMySelf() = 0;
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
