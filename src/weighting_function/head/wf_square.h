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



#ifndef __WF_SQUARE_H
#define __WF_SQUARE_H

#include "weighting_function.h"

/* --- Square --- */

class Square: public WeightingFunction
{
 protected:
  double m_factor;

  void init();

  virtual void setup();

 public:
  Square(Node *parent);
  virtual ~Square();

  virtual double interpolate(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    double d = m_cutoff;
  
    // since the function may be called with r=NULL, next line 
    // checks that; if it's so, we assume r->abs() = 0, which, 
    // e.g., for particle self-contributions is reasonable
    if(r) d -= r->abs();
    
    return m_factor * d * d;
  }

  virtual double weight(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
/*    throw gError
      ("Square::weightWithDist",
       "Please use a different weighting function. This one is not "
       "differential at r=rc.");

    return 0.;*/
    //MSG_INFO("Square::weight", "I AM HERE");
    if(!r) throw gError("Square::weight", "called with r=NULL !!! Contact the programmers.");
    
    double d = r->abs();
    if(d == 0.) throw gError("Square::weight", "called with r->abs() = 0 !!! Contact the programmers.");
    return 2*m_factor*(m_cutoff/d - 1);
  }

  virtual point_t localGradient(const Pairdist *r, const point_t &p, double *dist_from_wall = NULL) const {
    throw gError
      ("Square::weightWithDist",
       "Please use a different weighting function. This one is not "
       "differential at r=rc.");

    point_t g = {{{ 0, 0, 0 }}};
    return g;
  }
};

#endif
