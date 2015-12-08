/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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


#ifndef __TRIPLET_CALC_ANGULAR_FALPHA_H
#define __TRIPLET_CALC_ANGULAR_FALPHA_H

#include "triplet_calculator.h"

class Pairdist;
class ColourPair;

/*!
 * \a TripletCalculator for cached properties computed during a loop over bonded triplets
 * This \a Symbol specifically caches angular forces on particles derived 
 * from the potential 0.5*k_a*(a-a_0)^2 with force constant "k_a", 
 * "a" = PI - triplet angle and equilibrium angle "a_0" 
 */
class TripletCalcAngularFalpha : public TripletCalculator
{
protected:
  
  /*! force constant*/
  double m_k;
  /*! equilibrium angle*/
  double m_thetaEq;
  /*! helper: cosine of equilibrium angle*/
  double m_cosEq;
  /*! helper: Pi minus equilibrium angle*/
  double m_angleEq;
  
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   */
  virtual TripletCalcAngularFalpha* copyMySelf()
  {
    return new TripletCalcAngularFalpha(*this);
  }


public:
  /*!
 * Constructor
   */
  TripletCalcAngularFalpha(string symbol);
  
  /*!
   * Constructor
   */
  TripletCalcAngularFalpha(/*Node*/Simulation* parent);
   
  /*!
   * Destructor
   */
  virtual ~TripletCalcAngularFalpha() {
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
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
  bool findStage();

    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
  bool findStage_0();


};


#endif
