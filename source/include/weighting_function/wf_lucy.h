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



#ifndef __WF_LUCY
#define __WF_LUCY

#include "weighting_function.h"

/*!
 * Implementation of the Lucy weighting function
 * without walls
 */

class Lucy: public WeightingFunction
{
 protected:

  /*!
   * Precomputed factor for the interpolation function
   */
  double m_factor_i;

  /*!
   * Precomputed factor for the weighting function
   */
  double m_factor_w;

  /*!
   * Initialize property list
   */
  void init();

  /*!
   * Compute factors (i.e., \a m_factor_i and \a m_factor_w)
   */
  virtual void setup();

 public:

  /*!
   * Constructor.
   */
  Lucy(Node *parent);

  /*!
   * Destructor->
   */
  virtual ~Lucy();


  /*!
   * Compute interpolation function (i.e., Lucy function)
   * @param r pair distance
   * @param p position of the first particle (ignored)
   */
  virtual double interpolate(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    if (!r)
      return m_factor_i * pow(m_cutoff, 4);
    else
      return m_factor_i * (m_cutoff + 3*r->abs())*pow(m_cutoff - r->abs(), 3);
  }

  /*!
   * Compute weighting function (i.e., Derivative of Lucy function)
   * @param r pair distance
   * @param p position of the first particle (ignored)
   */
  virtual double weight(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    if (!r)
      return m_factor_w * pow(m_cutoff, 2);

    else
       return m_factor_w * pow(m_cutoff - r->abs(), 2);

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
