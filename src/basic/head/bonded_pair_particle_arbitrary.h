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


#ifndef __BONDED_PAIR_PARTICLE_ARBITRARY_H
#define __BONDED_PAIR_PARTICLE_ARBITRARY_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
#include "bonded_pair_particle_calc.h"
#include "colour_pair.h"
#include "function_pair.h"


/*!
 * Parent class with functions to compute completely user-defined properties for the
 * particles,  which need summation over bonded pairs. 
 */
class BondedPairParticleArbitrary : public BondedPairParticleCalc
{
 protected:
  
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
  
  /*!
   * Initialise the property list
   */
  virtual void init();

  /*!
   * The returned string contains those terms from runtime compiled expressions, 
   * which should be ignored when determining the stage. The expressions are separated by " | ".
   * An "empty" string must have the form "---".
   */
  virtual string usedSymbolsIgnoredForStaging() const {
    return m_oldSymbols;
  }
  
  /*!
   * Helper function for polymorphic copying
   */
  virtual ValCalculator* copyMySelf() = 0;
  
  /*!
   * copies the members of this class to \a vc
   */
  virtual void copyMembersTo(ValCalculator* vc);

  
 public:
  
  /*!
   * Constructor for the \a Node hierarchy
   */
  BondedPairParticleArbitrary(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~BondedPairParticleArbitrary();
  
#ifdef _OPENMP
  /*!
   * Merge the copies of all threads together
   */
  virtual void mergeCopies(ColourPair* cp, int thread_no) = 0;
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
