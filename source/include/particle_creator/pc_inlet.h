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



#ifndef __PC_INLET_H
#define __PC_INLET_H

#include "phase.h"
#include "pc_with_rng.h"
#include "boundary.h"
#include "manager_cell.h"
#include "particle_list.h"
#include "wall_triangle.h"


//---- Classes ----

class ParticleCreatorInlet: public ParticleCreatorWithRngPCalc
{
protected:
  /*!
  * particle density to be achieved in the inlet
  */
  double m_density;
  
  /*!
  * If the particle density is increased with a time ramp, 
  * this is the density we start from
  */
  double m_initDensity;
  
  /*!
  * number of timesteps for the time ramp
  */
  size_t m_nSteps;

  /*!
  * If the particle density is increased with a time ramp, this list contains the 
  * particles to be added per timestep. It will be popped until it is empty
  */
  list<size_t> m_createList;
  
  /*!
  * This is calculated from information given by a BoundaryWithInlet 
  */
  double m_inlet_volume, m_inlet_surface, m_inlet_length;
  point_t m_inlet_normal;

  /*! 
  * Store the face walls including their surface 
  */
  vector< pair<WallTriangle*, double> > m_walls;

  region_t *m_region;

  void init();

  virtual void scaleVels();

  virtual point_t randomPosition();

public:
  ParticleCreatorInlet(Boundary *boundary);

  virtual void setup();

  virtual void setRegion(region_t *r) {
    m_region = r;
  }

  virtual void createParticles();

  /*!
   * Create particles during the run
   */
  virtual void createMoreParticles()
  {
    if(!m_createList.empty())
    {
      size_t toCreate = m_createList.front();
      for(size_t i = 0; i < toCreate; ++i) createParticle(true);
      m_createList.pop_front();
    }
  }
  
  virtual void createParticle(bool assign_to_cell = true);
  virtual void deleteParticle();
};

#endif
