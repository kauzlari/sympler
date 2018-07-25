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

// #define VC_MAX_STAGE 1

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


#endif
