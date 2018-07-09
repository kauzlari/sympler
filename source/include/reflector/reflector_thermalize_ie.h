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



#ifndef __REFLECTOR_THERMALIZE_IE_H
#define __REFLECTOR_THERMALIZE_IE_H

#include "random.h"
#include "reflector_with_rng.h"
#include "function_fixed.h"
#include "reflector_thermalize.h"

/*!
 * The reflector reflects particles by drawing a new velocity that corresponds
 * to the temperature of the wall.
 *
 * Specifically, the velocity component perpendicular to the wall is drawn
 * from a Raleigh distribution and the velocity components parallel to the
 * wall are drawn from a normal distribution (for the corresponding wall temperature).
 * Additionally, the internal energy of the particle is set to the wall temperature
 * (without fluctuations).
 */
class ReflectorThermalizeInternalEnergy: public ReflectorWithRng
{
protected:
  /*!
   * Temperature of the wall (can be given in terms of the hit position)
   */
  FunctionFixed m_temperature;

  /*!
   * Square root of the temperature of the wall (can be given in terms of the hit position)
   */
  FunctionFixed m_temperature_sqrt;

  /*!
   * Internal energy to be set for reflected particles (can be given in terms of the hit position)
   */
  FunctionFixed m_internal_energy;

  /*!
   * The offset of the internal energy variable for each color
   */
  vector<size_t> m_ie_offset;
    
  /*!
   * Initialize the property list
   */
  void init();
    
public:
  /*!
   * Constructor
   * @param boundary The boundary this reflector is assigned to. Currently not used.
   */
  ReflectorThermalizeInternalEnergy(/*Wall *wall*/ Boundary* boundary);

  /*!
   * Destructor
   */
  virtual ~ReflectorThermalizeInternalEnergy();

  /*!
   * Check the temperature function and get the internal energy offsets
   */
  virtual void setup();
    
  /*!
   * Reflect a particle by choosing the reflection angle equal to the angle of incidence
   * @param p The particle that hit the wall
   * @param hit_pos Position, where the particle hit the wall
   * @param normal Normal vector to surface at the position where the particle hit the wall
   * @param in_plane A vector perpendicular to the surface at the position where the particle hit the wall
   */
  virtual void reflect
      (Particle *p, point_t& r, point_t& v, const point_t &hit_pos, const point_t &normal, const point_t &in_plane) {
    double vel_perp, vel_para1, vel_para2, ts;
    point_t in_plane2;

    ts = m_temperature_sqrt(hit_pos.x, hit_pos.y, hit_pos.z);
		
    in_plane2 = normal.cross(in_plane);
		
    //		vel_perp = gsl_ran_rayleigh(s_rng, m_temperature_sqrt);
    //		vel_para1 = s_rng_normal.nrandf()*m_temperature_sqrt;
    //		vel_para2 = s_rng_normal.nrandf()*m_temperature_sqrt;
    vel_perp = m_rng.rayleigh(ts);
    vel_para1 = m_rng.normal(1)*ts;
    vel_para2 = m_rng.normal(1)*ts;
		
    /*p->*/r = hit_pos + c_rt_disp_eps*normal;
    /*p->*/v = vel_perp*normal + vel_para1*in_plane + vel_para2*in_plane2;
    p->tag.doubleByOffset(m_ie_offset[p->c]) = m_internal_energy(hit_pos.x, hit_pos.y, hit_pos.z);

    ++s_n_hits;
  }
};

#undef EPSILON

#endif
