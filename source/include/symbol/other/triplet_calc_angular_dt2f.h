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


#ifndef __TRIPLET_CALC_ANGULAR_DT2F_H
#define __TRIPLET_CALC_ANGULAR_DT2F_H

#include "triplet_calculator.h"

class Pairdist;
class ColourPair;

/*!
 * \a TripletCalculator for cached properties computed during a loop over bonded triplets
 * This \a Symbol specifically caches the 0th, 1st, and 2nd time derivatives of angular forces on particles derived 
 * from the potential 0.5*k_a*(cosa-cosa_0)^2 with force constant "k_a", 
 * cosine of [PI - triplet angle] "cosa" and cosine of equilibrium angle "cosa_0" 
 */
class TripletCalcAngularDt2F : public TripletCalculator
{
protected:

  /*! Should we assume 2D? Influences only the usage of the usage of the function g_idMinusEqTensor2D instead of the 3D version.*/
  bool m_2D;  

  /*! force constant*/
  double m_k;
  /*! equilibrium angle*/
  double m_thetaEq;
  /*! helper: cosine of equilibrium angle*/
  double m_cosEq;

  ////////// FORCES NOT DONE HERE, JUST THE TWO DERIVATIVES!!! ///////////
  // So we need one additional symbol
#if 0
  /*!
   * Memory offsets to the symbol in the \a Particle 's tag this calculator uses to store the computed FIRST DERIVATIVE of the force, one for each of the possibly different species 
   */
  size_t m_slotsDtF[3];
  
  /*!
   * Name of the symbol storing the first derivative of the force
   */
  string m_DtFname;
#endif

  /*!
   * Memory offsets to the symbol in the \a Particle 's tag this calculator uses to store the computed SECOND DERIVATIVE of the force, one for each of the possibly different species 
   */
  size_t m_slotsDt2F[3];

  /*!
   * Name of the symbol storing the second derivative of the force
   */
  string m_Dt2Fname;
  
  /*!
   * Memory offsets to the symbol in the \a Particle 's tag this calculator uses to read the force needed to compute its 2nd derivative; one for each of the possibly different species 
   */
  size_t m_slotsFin[3];

  /*!
   * Name of the symbol storing the second derivative of the force
   */
  string m_FinName;
  
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   */
  virtual TripletCalcAngularDt2F* copyMySelf()
  {
    return new TripletCalcAngularDt2F(*this);
  }

  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
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
    usedSymbols.push_back(m_FinName);
  }

    
public:
  /*!
 * Constructor
   */
  TripletCalcAngularDt2F(string symbol);
  
  /*!
   * Constructor
   */
  TripletCalcAngularDt2F(/*Node*/Simulation* parent);
   
  /*!
   * Destructor
   */
  virtual ~TripletCalcAngularDt2F() {
  }

  /*!
   * Compute cached properties for \a triplet_t \a tr
   * @param tr The \a triplet_t for which to compute the cached properties
   */
#ifndef _OPENMP
  virtual void compute(triplet_t* tr);
#else
  virtual void compute(triplet_t* tr/*, size_t thread_no*/ /*FIXME: parallelise!*/);
#endif


#ifdef _OPENMP
  /*!
   * Merging the copies among different threads (processors) together
   */
    virtual void mergeCopies(size_t thread_no);

    virtual int setNumOfDoubles();

#endif

    /*!
     * setup this calculator
     */
    virtual void setup();
    
    /*!
     * run those setups requiring that \a Particle s have been created
     */
    virtual void setupAfterParticleCreation();
    
    /*!
     * Return the name of the computed symbols to be used in other expressions.
     */
    virtual list<string> mySymbolNames()
    {
      list<string> temp;
      assert(temp.empty());
      temp.push_back(m_symbolName);
      temp.push_back(m_Dt2Fname);
      return temp;
    }
    
};

#endif
