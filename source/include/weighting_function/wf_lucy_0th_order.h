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



#ifndef __WF_LUCY_0TH_ORDER
#define __WF_LICY_0TH_ORDER

#include "weighting_function.h"


/*!
 * Implementation of the Lucy weighting function
 * with walls
 */
class LucyWithWall0thOrder: public WeightingFunctionWithWall
{
 protected:
  /*!
   * 1/cutoff^7
   */
  double m_rc7;
  
  /*!
   * Cached factor
   */
  double m_factor;

  /*!
   * Compute the first order correction to the weighting function
   * if close to a wall
   * @param zp Distance to the next wall normalized to the cutoff radius
   */
  double A(double zp) const {
    assert(zp >= 0);

    return 5./
      (M_PI*m_rc7*
       (8./21 + zp - (5./3)*pow(zp,3) +
	3*pow(zp,5) - (8./3)*pow(zp,6) + (5./7)*pow(zp,7)
	)
       );
  }

  /*!
   * Compute the negative of the derivative of first order
   * correction to the weighting function
   * if close to a wall
   * @param zp Distance to the next wall normalized to the cutoff radius
   */
  double neg_dAdz(double zp) const {
    assert(zp >= 0);
    
    return 5.*
      (1-5*pow(zp,2)+15*pow(zp,4)-16*pow(zp,5)+5*pow(zp,6))	 
      /
      (M_PI*m_rc7*m_cutoff*
       pow(8./21 + zp - (5./3)*pow(zp,3) + 
	   3*pow(zp,5) - (8./3)*pow(zp,6) + (5./7)*pow(zp,7), 2));
  }

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
   * Constructor->
   */
  LucyWithWall0thOrder(Node *parent);

  /*!
   * Destructor->
   */
  virtual ~LucyWithWall0thOrder();

  /*!
   * Compute interpolation function (i.e., Lucy function)
   * @param r Pair distance
   * @param normal Normal vector to the surface
   * @param dist Distance to the next wall normalized to the cutoff radius
   */
  virtual double interpolateWithDist(const Pairdist *r, const point_t &normal, double dist) const {
    if (!r)
      return A(dist)*pow(m_cutoff, 4);
    else
      return A(dist)*(m_cutoff + 3*r->abs())*pow(m_cutoff - r->abs(), 3);
  }

  /*!
   * Compute weighting function (i.e., Derivative of Lucy function)
   * @param r Pair distance
   * @param normal Normal vector to the surface
   * @param dist Distance to the next wall normalized to the cutoff radius
   */
  virtual double weightWithDist(const Pairdist *r, const point_t &normal, double dist) const {
    if (!r)
      return 12*A(dist)*pow(m_cutoff, 2);
    else
      return 12*A(dist)*pow(m_cutoff - r->abs(), 2);
  }

  /*!
   * Compute local gradient functionl
   * @param r Pair distance
   * @param normal Normal vector to the surface
   * @param dist Distance to the next wall normalized to the cutoff radius
   */
  virtual point_t localGradientWithDist(const Pairdist *r, const point_t &normal, double dist) const {
    if (!r)
      return normal*neg_dAdz(dist)*pow(m_cutoff, 4);
    else
      return normal*neg_dAdz(dist)*(m_cutoff + 3*r->abs())*pow(m_cutoff - r->abs(), 3);
  }

  /*!
   * For a constant surface entropy factor return the gradient
   * grad(\int d^2 r_s W(r_s-r_i, r_i))
   * @param normal The normal vector to the wall
   * @param dist The distance from the wall
   */
  virtual point_t gradientSurfaceWeightWithDist(const point_t &normal, double dist) const {
    point_t g = {{{ 0, 0, 0 }}};

    if (dist < 1) {
      g =
	-
	(M_PI*pow(m_cutoff, 6)/5. * neg_dAdz(dist)*pow(1-dist, 4)*(1+4*dist+5*dist*dist)
	 +
	 2*M_PI*pow(m_cutoff, 5) * A(dist) * dist*pow(1-dist, 3)*(1+3*dist)
	 )
	* normal;
    }

    return g;
  }
};

#endif
