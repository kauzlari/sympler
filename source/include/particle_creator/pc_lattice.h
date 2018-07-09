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


#ifndef __PC_LATTICE_H
#define __PC_LATTICE_H

#include "phase.h"
#include "pc_with_rng.h"
#include "boundary.h"
#include "particle_list.h"


//---- Classes ----

class ParticleCreatorLattice: public ParticleCreatorWithRngF
{
protected:
  /*!
  *Number of lattice points in each direction 
  */
  int_point_t m_nlattice_points;
  
  /*!
  * nearest neighbour distance of the particles
  */
  double m_distance;
  /*!
  * Number density of the particles
  */
  double m_density/*, m_temperature*/;

//   /* If randomize = false then seed = 1 (for rand) */
//   bool m_randomize;

  /*! If force_boundary_size is on, the boundary will be resized to
     * fit the specified information.
  */
  bool m_force_boundary_size;

  /*!
  * Helpers
  */
  bool m_lattice_defined, m_dd_defined;
      
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
  * Scale the particle velocities according to \a m_temperature
  */
  virtual void scaleVels();
  
  virtual void computeDistance();

  virtual double spacing(const double& distance) const;

  public:
  /*!
  * Constructor
  */
  ParticleCreatorLattice(Boundary *boundary);
	
  /*!
  * Check whether the defined bounding box must be resised by 
  * this \a ParticleCreator
  */
	virtual void adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend);

 /*!
 * Create the particles
 */
 virtual void createParticles();
		
};

#endif
