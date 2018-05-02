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



#ifndef __WF_LUCY_1ST_ORDER
#define __WF_LICY_1ST_ORDER

#include "weighting_function.h"


/*!
 * Implementation of the Lucy weighting function
 * with walls
 */

class LucyWithWall1stOrder: public WeightingFunctionWithWall
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
   * Compute the 0th order correction to the weighting function
   * if close to a wall
   * @param dist Pair distance
   */
  double A0(double dist) const {
    double z0 = fabs(dist);

    /* Don't ask....
       If you do, it has been autocomputed using Mathematica
    */
#define Power(x, y) pow(x,y)
    return 
      (6720*
       (2 + Power(z0,3)*
	(21 + Power(z0,2)*
	 (-63 + Power(z0,2)*
	  (135 + 7*z0*
	   (-18 + 5*z0)
	   )
	  )
	 )
	)
       )
      /
      (M_PI*m_rc7*
       (-27*Power(-1 + z0,10)*Power(5 + z0*(25 + z0*(47 + 35*z0)),2) + 
	64
	*
	(8 + 21*z0 - 35*Power(z0,3) + 63*Power(z0,5) - 56*Power(z0,6) + 15*Power(z0,7))
	*
	(2 + Power(z0,3)*
	 (21 + Power(z0,2)*
	  (-63 + Power(z0,2)*
	   (135 + 7*z0*(-18 + 5*z0))
	   )
	  )
	 )
	)
       );
#undef Power
  }


  /*!
   * Compute the 1st order correction to the weighting function
   * if close to a wall
   * @param dist Pair distance
   */
  double A1(double dist) const {
    double z0 = fabs(dist);

    /* Don't ask....
       If you do, it has been autocomputed using Mathematica
    */
#define Power(x, y) pow(x,y)
    return
      (7560*Power(-1 + z0,5)*
       (5 + z0*
	(25 + z0*
	 (47 + 35*z0)
	 )
	)
       )
      /
      (M_PI*m_rc7*m_cutoff*
       (-27*Power(-1 + z0,10)*Power(5 + z0*(25 + z0*(47 + 35*z0)),2) + 
	64*
	(8 + 21*z0 - 35*Power(z0,3) + 63*Power(z0,5) - 56*Power(z0,6) + 15*Power(z0,7))
	*
	(2 + Power(z0,3)
	 *
	 (21 + Power(z0,2)*
	  (-63 + Power(z0,2)*
	   (135 + 7*z0*
	    (-18 + 5*z0)
	    )
	   )
	  )
	 )
	)
       );
#undef Power
  }


  /*!
   * Compute the negative of the derivative of 0th order
   * correction to the weighting function
   * if close to a wall
   * @param z0 Distance to the next wall normalized to the cutoff radius
   */
  double dA0dz(double z0) const {
    assert(z0 >= 0);
    
#define Power(x, y) pow(x,y)
    return 
(-141120*Power(-1 + z0,4)*Power(16 + 45*z0 - 84*Power(z0,3) + 126*Power(z0,5) - 
       180*Power(z0,7) + 144*Power(z0,8) - 35*Power(z0,9),2)*(1 + z0*(4 + 5*z0)))/
   (M_PI*m_rc7*m_cutoff*Power(349 + z0*(2688 + z0*(7560 + 
             z0*(6272 + z0*(-11844 + 
                   z0*(-24192 + z0*
                       (4760 + z0*(36480 + 
                            z0*(8190 + 
                              z0*(-33152 + 
                              z0*(-12936 + 
                              z0*(24192 + 
                              z0*(8540 + 
                              z0*(-24192 + z0*(15336 + 35*z0*(-128 + 15*z0)))))))))))
			    )))),2));
#undef Power
  }


  /*!
   * Compute the negative of the derivative of 1st order
   * correction to the weighting function
   * if close to a wall
   * @param z0 Distance to the next wall normalized to the cutoff radius
   */
  double dA1dz(double z0) const {
    assert(z0 >= 0);
    
#define Power(x, y) pow(x,y)
    return
(-423360*Power(-1 + z0,4)*(-16 - 45*z0 + 84*Power(z0,3) - 126*Power(z0,5) + 
       180*Power(z0,7) - 144*Power(z0,8) + 35*Power(z0,9))*(1 + z0*(4 + 5*z0))*
     (15 + z0*(64 + 84*z0 - 70*Power(z0,3) + 84*Power(z0,5) - 64*Power(z0,6) + 
          15*Power(z0,7))))/
   (M_PI*m_rc7*m_cutoff*m_cutoff*Power(349 + z0*(2688 + z0*(7560 + 
             z0*(6272 + z0*(-11844 + 
                   z0*(-24192 + z0*
                       (4760 + z0*(36480 + 
                            z0*(8190 + 
                              z0*(-33152 + 
                              z0*(-12936 + 
                              z0*(24192 + 
                              z0*(8540 + 
                              z0*(-24192 + z0*(15336 + 35*z0*(-128 + 15*z0)))))))))))
			    )))),2));
#undef Power
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
  LucyWithWall1stOrder(Node *parent);

  /*!
   * Destructor->
   */
  virtual ~LucyWithWall1stOrder();

  /*!
   * Compute interpolation function (i.e., Lucy function)
   * @param r pair distance
   * @param normal Normal vector to the surface
   * @param dist distance of the first particle from the wall
   */
  virtual double interpolateWithDist(const Pairdist *r, const point_t &normal, double dist) const {
    if (!r)
      return A0(dist)*pow(m_cutoff, 4);
    else
      return
	(A0(dist)+A1(dist)*(r->cartesian()*normal))
	*
	(m_cutoff + 3*r->abs())*pow(m_cutoff - r->abs(), 3);
  }

  /*!
   * Compute weighting function (i.e., Derivative of Lucy function)
   * @param r pair distance
   * @param normal Normal vector to the surface
   * @param dist distance of the first particle from the wall
   */
  virtual double weightWithDist(const Pairdist *r, const point_t &normal, double dist) const {
    /*    if (!r)
      return
	-
	A1(dist)*pow(m_cutoff, 4)/r->abs()
	+
	12*A1(dist)*(r->cartesian()*normal)*pow(m_cutoff - r->abs(), 2);
	else {*/
    assert(r);
    assert(r->abs() > 0);

    double a0 = A0(dist);
    double a1 = A1(dist);

    return
      -
      (
       a1*(m_cutoff + 3*r->abs())*pow(m_cutoff - r->abs(), 3)/r->abs()
       -
       12*(a0+a1*(r->cartesian()*normal))*pow(m_cutoff - r->abs(), 2)
      );
  }

  /*!
   * Compute local gradient function
   * @param r pair distance
   * @param normal Normal vector to the surface
   * @param dist distance of the first particle from the wall
   */
  virtual point_t localGradientWithDist(const Pairdist *r, const point_t &normal, double dist) const {
    if (!r)
      return
	-
	dA0dz(dist)*normal
	*
	pow(m_cutoff, 4);
    else
      return
	-
	(dA0dz(dist) + dA1dz(dist)*(r->cartesian()*normal))*normal
	*
	(m_cutoff + 3*r->abs())*pow(m_cutoff - r->abs(), 3);
  }
};

#endif
