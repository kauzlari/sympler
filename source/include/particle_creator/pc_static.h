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


/* $ */

#ifndef __PC_STATIC_H
#define __PC_STATIC_H

#include "phase.h"
#include "pc_with_rng.h"
#include "boundary.h"
#include "particle_list.h"


/*! 
 *---- Classes ----
*/

class ParticleCreatorStatic: public ParticleCreatorWithRngPCalc
{
protected:

  /*!
 * Should more particles be created DURING the simulation?
  */
  bool m_more;
  /*!
   * Should the particle creator explicitely check whether the 
   * particles that should be created are really inside?
   */
  bool m_checkInside;
  
  /*!
 * Time between creation of more particles. Will be ignored if m_more = false.
  */
  double m_dtForMore;

  /*!
 * Time. The first particles are created at the time t = timeOffset. Will be ignored if m_more = false.
  */
  double m_timeOffset;

  /*!
   * Internal variable next absolute time for creation of more particles. Will be ignored if m_more = false.
   */
  double m_nextTimeForMore;
    
  /*!
   * Should we create the inscribed ellipsoid?
   */
  bool m_ellipsoid;

  /*!
 * The points (particles) that are defined in every direction
 */

  int_point_t m_nstatic_points;

    /*!
   * Position of the corners
   */
  point_t m_corners[3];

  /*!
   * Vertices of the corners.
   */
  int m_corners_v[3];
  
//   point_t m_corner1, m_corner2/*, point1, point2*/;
  
/*! 
 *The cuboid that defines the staticBox  
*/  
 
  cuboid_t m_static_box;
    
  double m_distance, m_density/*, m_temperature*/;

//   bool m_randomize;
 
  bool m_force_boundary_size;

  point_t m_offset, m_size;
  
  bool m_static_defined, m_ddstatic_defined;
    
  void init();

  virtual void scaleVels();

public:


  ParticleCreatorStatic(Boundary *boundary);
  
  virtual const cuboid_t &staticBox() const {
    return m_static_box;
  }
  
  /*! 
   * Creates a static Box and takes care that it's created in the boundary
   */

  virtual void adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend);

		
  virtual void createParticles();

  /*!
   * Create particles during the run
   */
  virtual void createMoreParticles();

};

#endif
