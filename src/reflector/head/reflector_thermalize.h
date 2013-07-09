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



#ifndef __REFLECTOR_THERMALIZE_H
#define __REFLECTOR_THERMALIZE_H

#include "random.h"
#include "reflector_with_rng.h"
#include "boundary.h"
#include "phase.h"
#include "simulation.h"
#include "function_fixed.h"

extern double c_rt_disp_eps;

/*!
 * The reflector reflects particles by drawing a new velocity that corresponds
 * to the temperature of the wall.
 *
 * Specifically, the velocity component perpendicular to the wall is drawn
 * from a Raleigh distribution and the velocity components parallel to the
 * wall are drawn from a normal distribution (for the corresponding wall temperature).
 */
class ReflectorThermalize: public ReflectorWithRng
{
 protected:
  /*!
   * Temperature of the wall
   */
  double m_temperature;

  /*!
   * Square root of the temperature of the wall
   */
  double m_temperature_sqrt;

  /*!
   * A bias velocity by which the velocity component parallel to the wall
   * is shifted.
   */
//   double m_bias_velocity;
  FunctionFixed m_bias_velocity;
  
  /*!
   * Determines the orientation of the bias velocity in dependence of 
   * \a m_rotational_bias
   */
  point_t m_bias_orientation;
    

  /*!
   * Boolean switch, needed to determine, the applied bias velocity. 
   * False-case: \a m_bias_velocity * \a m_bias_orientation
   * True-case: \a m_bias_velocity * (normal x \a m_bias_orientation), where 
   * normal is the surface normal of the wall.
   */
  bool m_rotational_bias;

  /*!
  * If this attribute is set to 'true', the tracking of the particles 
  * trajectory will be aborted after the first hit of the particle at a wall
  */
  bool m_oneHit;
    
  /*!
   * Initialize the property list
   */
  void init();
    
 public:
  /*!
   * Constructor
   * @param boundary The boundary this reflector is assigned to. Currently not used.
   */
   ReflectorThermalize(/*Wall *wall*/ Boundary* boundary);

  /*!
   * Constructor
   * @param boundary The boundary this reflector is assigned to. Currently not used.
   * @param temperature The temperature of the wall
   */
   ReflectorThermalize(/*Wall *wall*/ Boundary* boundary, double temperature);

  /*!
   * Destructor
   */
  virtual ~ReflectorThermalize();

  /*!
   * Initialize the cached square root of the temperature
   * and normalize the second normal.
   */
  virtual void setup();
    
  /*!
   * Reflect a particle by choosing a new velocity
   * @param p The particle that hit the wall
   * @param hit_pos Position, where the particle hit the wall
   * @param normal Normal vector to surface at the position where the particle hit the wall
   * @param in_plane A vector perpendicular to the surface at the position where the particle hit the wall
   */
  virtual void reflect
      (Particle *p, point_t& r, point_t& v, const point_t &hit_pos, const point_t &normal,
     const point_t &in_plane) {

    //    MSG_DEBUG("ReflectorThermalize::reflect", "hit_pos=(" + ObjToString(hit_pos.x) + ", " + ObjToString(hit_pos.y) + ", " + ObjToString(hit_pos.z) + ")");
    
    double vel_perp, vel_para1, vel_para2;
    point_t in_plane2;
    point_t bias = {{{0, 0, 0}}};

//     Boundary* b = (Boundary*) m_parent;
//     Phase* p = (Phase*) (b->parent());
//     Simulation* s = (Simulation*) (p->parent());
        double time = 
          (
            (Simulation*)
            (
              (
                (Phase*)
                (
                  (
                    (Boundary*) m_parent
                  )->parent()
                )
              )->parent()
            )
          )->controller()->time();
          
    if(!m_rotational_bias)
      bias = m_bias_velocity(time, hit_pos.x, hit_pos.y, hit_pos.z)*m_bias_orientation;
    else
    
      bias = m_bias_velocity(time, hit_pos.x, hit_pos.y, hit_pos.z)*normal.cross(m_bias_orientation);

    //    MSG_DEBUG("ReflectorThermalize::reflect", "bias = " << bias << ", name = " << m_reflector_name);

    in_plane2 = normal.cross(in_plane);
		
    vel_perp = m_rng.rayleigh(m_temperature_sqrt);
    vel_para1 = m_rng.normal(1)*m_temperature_sqrt;
    vel_para2 = m_rng.normal(1)*m_temperature_sqrt;
		
    r = hit_pos + c_rt_disp_eps*normal;
    v = vel_perp*normal + vel_para1*in_plane + vel_para2*in_plane2 + bias;
    
    if(m_oneHit) p->dt = 0;
    
    ++s_n_hits;
  }
};

#undef EPSILON

#endif
