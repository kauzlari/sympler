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



#ifndef __WF_INPUTWF
#define __WF_INPUTWF

#include "weighting_function.h"
#include "function_fixed.h"

/*!
 * Implementation of the weighting function definable by input file
 * without walls
 */

class InputWF: public WeightingFunction
{
 protected:

  /*!
   * Entered interpolation function W(r) 
   */
  FunctionFixed m_weightFunct;

  /*!
   * Entered weighting function derivative F(r) = -W'(r)/r
   */
  FunctionFixed m_Deriv_wF;

  /*!
   * Self contribution factor for interpolation function
   */
  FunctionFixed m_selfContrib;

//  /*!
//   * Self contribution factor for weighting function
  //   */
//  FunctionFixed m_Deriv_sC;


  /*!
   * Initialize property list
   */
  void init();

  /*!
   * Get functions (i.e., m_weightFunct, m_Deriv_wF, m_selfContrib & m_Deriv_sC)
   */
  virtual void setup();

 public:

  /*!
   * Constructor.
   */
  InputWF(Node *parent);

  /*!
   * Destructor->
   */
  virtual ~InputWF();


  /*!
   * Compute interpolation function (i.e., input function)
   * @param r pair distance
   * @param p position of the first particle (ignored)
   */
  virtual double interpolate(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    if (!r)
      return m_selfContrib(m_cutoff);

    else
      return m_selfContrib(m_cutoff) * m_weightFunct(r->abs());
  }

  /*!
   * Compute weighting function (i.e., Derivative of input function)
   * @param r pair distance
   * @param p position of the first particle (ignored)
   */
  virtual double weight(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    if 
      (!r)throw gError("InputWF::weight", "called with r=NULL !!! Contact the programmers.");
     
    else
      return /*m_Deriv_sC(m_cutoff) * */m_Deriv_wF(r->abs());
  }

  /*!
   * Compute local gradient functionl
   * @param r pair distance
   * @param p position of the first particle (ignored)
   */
  virtual point_t localGradient(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    point_t g = {{{ 0, 0, 0 }}};

    return g;
  }
};

#endif
