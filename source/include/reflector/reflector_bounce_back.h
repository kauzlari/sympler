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


#ifndef __REFLECTOR_BOUNCE_BACK_H
#define __REFLECTOR_BOUNCE_BACK_H 

#include "random.h"
#include "reflector_with_rng.h"

extern double c_rm_disp_eps;

/*!
 * This reflector reflects particles by inverting their velocity
 */
class ReflectorBounceBack: public Reflector
{
protected:
  /*!
   * Initialize the property list
   */
  void init();
    
public:
  /*!
   * Constructor
   * @param boundary The boundary this reflector is assigned to. Currently not used.
   */
  ReflectorBounceBack(/*Wall *wall*/ Boundary* boundary);

  /*!
   * Destructor
   */
  virtual ~ReflectorBounceBack();
    
  /*!
   * Reflect a particle by inverting the velocity
   * @param p The particle that hit the wall
   * @param hit_pos Position, where the particle hit the wall
   * @param normal Normal vector to surface at the position where the particle hit the wall
   * @param in_plane A vector perpendicular to the surface at the position where the particle hit the wall
   */
  virtual void reflect
      (Particle *p, point_t& r, point_t& v, const point_t &hit_pos,
     const point_t &normal, const point_t &in_plane) {
		
    /*p->*/r = hit_pos + c_rm_disp_eps*normal;
    /*p->*/v = -v;

    ++s_n_hits;
  }
};

#endif
