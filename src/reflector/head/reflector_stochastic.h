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



#ifndef __REFLECTOR_STOCHASTIC_H
#define __REFLECTOR_STOCHASTIC_H

#include "random.h"
#include "reflector_with_rng.h"

extern double c_rs_disp_eps;

/*!
 * The reflector reflects particles by keeping their velocity
 * and choosing a random reflection angle.
 *
 * Specifically, the reflection angles are choosen from a uniform distribution.
 */
class ReflectorStochastic: public ReflectorWithRng
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
  ReflectorStochastic(/*Wall *wall*/ Boundary* boundary);
  
  /*!
   * Destructor
   */
  virtual ~ReflectorStochastic();
  
  /*!
   * Reflect a particle by choosing a random reflection angle
   * @param p The particle that hit the wall
   * @param hit_pos Position, where the particle hit the wall
   * @param normal Normal vector to surface at the position where the particle hit the wall
   * @param in_plane A vector perpendicular to the surface at the position where the particle hit the wall
   */
  virtual void reflect(Particle *p, point_t& r, point_t& v, const point_t &hit_pos, const point_t &normal, const point_t &in_plane) {
    double vel, vel_sin_phi, vel_perp, vel_para1, vel_para2, phi, theta;
    point_t in_plane2;
		
    in_plane2 = normal.cross(in_plane);

    phi = m_rng.uniform()*M_PI/2;
    theta = m_rng.uniform()*2*M_PI;

    vel = /*p->*/v.abs();
    vel_sin_phi = vel*sin(phi);

    vel_perp = vel*cos(phi);
    vel_para1 = vel_sin_phi*cos(theta);
    vel_para2 = vel_sin_phi*sin(theta);
		
    /*p->*/r = hit_pos + c_rs_disp_eps*normal;
    /*p->*/v = vel_perp*normal + vel_para1*in_plane + vel_para2*in_plane2;

    ++s_n_hits;
  }
};

#endif
