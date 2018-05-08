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



#ifndef __CBL_PAIR_PARTICLE_ARBITRARY_H
#define __CBL_PAIR_PARTICLE_ARBITRARY_H

#include "callable.h"

#include "function_pair.h"

/*!
 * \A Callable for cached properties computed during a loop over pairs.
 * Any subclass of abstract class \a CblPairParticleArbitrary computes 
 * a completely user-defined property for each particle, which needs 
 * pair summation.
 */
class CblPairParticleArbitrary : public Callable
{

 protected:

  /*!
   * The species of the ColourPair this \A Callable should belong to
   */
  pair<string, string> m_species;

  /*!
   * The tag offset of the data stored in a \a Particle
   */
  pair<size_t, size_t> m_slots;

  /*!
   * Cut-off radius for the pair summation
   */
  double m_cutoff;
  
  /*!
   * The name of the computed symbol to be used in other expressions
   */
  string m_symbolName;

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
   * Should all colour combinations be considered. This disables \a m_species.
   */
  bool m_allPairs;

  /*!
   * How does the pair expression behave under index interchange?
   * 1 means symmetric behaviour
   * -1 means antisymmetric behaviour
   */
  int m_symmetry;

  /*!
   * Is this \a Callable allowed to overwrite already existing symbols 
   * with name \a m_symbolName ?
   */
  bool m_overwrite;

  /*!
   * The datatype of this \a Callable (usually \a DataFormat::DOUBLE, 
   * POINT, TENSOR)
   */
  DataFormat::datatype_t m_datatype;

  /*!
   * The pair of colors of the \a Particle s to compute on
   */
  ColourPair *m_cp;
  
  /*!
   * Initialise the property list
   */
  virtual void init();

  
public:

  /*!
   * Constructor
   */
  CblPairParticleArbitrary(string symbol);

  /*!
   * Constructor
   */
  CblPairParticleArbitrary(/*Node*/Simulation* parent);

  /*!
   * Destructor
   */
  virtual ~CblPairParticleArbitrary() {
  }

  /*!
   * Setup this \a Callable
   */
  virtual void setup();

  /*!
   * Return the pair of species this \a Callable is computing for
   */
  const pair<string, string>& returnSpeciesPair() const {
    return m_species;
  }

  /*!
   * Return the name of the computed symbol
   */
  const string& returnSymbolName() const {
    return m_symbolName;
  }
  
  /*!
   * Return the specified symmetry of \a m_expression
   */
  size_t returnSymmetry() const {
    return m_symmetry;
  }
  
  /*!
   * Return \a m_cutoff
   */
  const double& returnCutoff() const {
    return m_cutoff;
  }
  
  /*!
   * Return the \a m_overwrite flag
   */
  bool returnOverwriteFlag() const {
    return m_overwrite;
  }

  /*!
   * Return the runtime compiled pair expression string
   */
  const string& returnPairExpression() const {
    return m_expression;
  }

  /*!
   * Return the runtime compiled expression string for 
   * \a m_1stparticleFactor
   */
  const string& return1stParticleExpression() const {
    return m_1stPExpression;
  }

  /*!
   * Return the runtime compiled expression string for 
   * \a m_2ndparticleFactor
   */
  const string& return2ndParticleExpression() const {
    return m_2ndPExpression;
  }

};

#endif

