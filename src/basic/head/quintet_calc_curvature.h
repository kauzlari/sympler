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


#ifndef __QUINTET_CALC_CURVATURE_H
#define __QUINTET_CALC_CURVATURE_H

#include "math.h"
#include "quintet_calculator.h"

class Pairdist;
class ColourPair;

/*!
 * \a TripletCalculator for cached properties computed during a loop over bonded triplets
 * This \a Symbol specifically caches angular forces on particles derived 
 * from the potential 0.5*k_a*(cosa-cosa_0)^2 with force constant "k_a", 
 * cosine of [PI - triplet angle] "cosa" and cosine of equilibrium angle "cosa_0" 
 */
class QuintetCalcCurvature : public QuintetCalculator
{
protected:
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   */
  virtual QuintetCalcCurvature* copyMySelf()
  {
    return new QuintetCalcCurvature(*this);
  }


public:
  /*!
 * Constructor
   */
  QuintetCalcCurvature(string symbol);
  
  /*!
   * Constructor
   */
  QuintetCalcCurvature(/*Node*/Simulation* parent);
   
  /*!
   * Destructor
   */
  virtual ~QuintetCalcCurvature() {
  }

  /*!
   * Compute cached properties for \a quintet_t \a q
   * @param q The \a quintet_t for which to compute the cached properties
   */
#ifndef _OPENMP
  virtual void compute(quintet_t* q);
#else
  virtual void compute(quintet_t* q/*, size_t thread_no*/ /*FIXME: parallelise!*/);
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
