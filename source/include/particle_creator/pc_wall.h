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


#ifndef __PC_WALL_H
#define __PC_WALL_H

#include "phase.h"
#include "boundary.h"
#include "particle_creator.h"
#include "function_fixed.h"

/* ---- ParticleCreatorWall ---- */

/*!
 * Generates particles behind walls, i.e., outside the domain,
 * a freely moving particle can reach.
*/
class ParticleCreatorWall: public ParticleCreator
{
protected:
  
  /*!
  * Sets the x-component of the position to the specified algebraic expression. See init() for the known symbols.
  */
  FunctionFixed m_posX;
  
    /*!
   * Sets the y-component of the position to the specified algebraic expression. See init() for the known symbols.
     */
  FunctionFixed m_posY;
  
    /*!
   * Sets the z-component of the position to the specified algebraic expression. See init() for the known symbols.
     */
  FunctionFixed m_posZ;
  
  /*!
  * Sets the x-component of the velocity to the specified algebraic expression. See init() for the known symbols.
  */
  FunctionFixed m_velX;
  
    /*!
   * Sets the y-component of the velocity to the specified algebraic expression. See init() for the known symbols.
     */
  FunctionFixed m_velY;
  
    /*!
   * Sets the z-component of the velocity to the specified algebraic expression. See init() for the known symbols.
     */
  FunctionFixed m_velZ;
  
  /*! 
  * Number of lattice points in each direction
  */
  math_vector_t<int> m_nlattice_points;

  /*!
  * Spacing of the lattice points.
  */
  double m_distance;
  
  /*!
  * Number density of the particles.
  */
  double m_density;

  /*!
   * Takes the cutoff from SImulation before the "skin" from the Verlet List is added
   */
  double myCutoff;
  
  /*!
   * Generate frozen particles at the 'upper end' of the x/y/z-direction in any case
  */
  bool_point_t m_forceMax;
  
  /*!
   * Generate frozen particles at the 'lower end' of the x/y/z-direction in any case
   */
  bool_point_t m_forceMin;
  
  /*!
   * Should the particles be frozen particles, i.e., all quantities fixed to
  * the initial value and unchangeable?
  */
  bool m_freeze;
		
  /*!
  * There may be only one ParticleCreatorWall in the simulation. This static
  * variable controls that. 
  * FIXME: eventually one PCWall PER SPECIES could become necessary!
  */
  static bool already;
     
  /*!
  * Initialisation of the \a PropertyList
  */
  void init();

  /*!
  * Setup for the \a ParticleCreatorWall
  */
  virtual void setup();

public:
  /*!
  * Constructor
  * @param boundary The \a Boundary, this \a ParticleCreator belongs to 
  */
  ParticleCreatorWall(Boundary *boundary);
		
  /*!
  * Destructor
  */
  virtual ~ParticleCreatorWall();
	
  /*!
  * In principle for adjusting the box size proposed by the boundary. But in fact, 
  * this ParticleCreator adjusts the box size only indirectly via the second and 
  * third argument
  * @param size the box size proposed by the \a Boundary , this 
  * \a ParticleCreator belongs to 
   * @param frameRCfront Should the 'lower end' of the simulation box be extended 
   * in order to put wall \a Particle s there?
   * @param frameRCend Should the 'upper end' of the simulation box be extended 
   * in order to put wall \a Particle s there?
  */
  virtual void adjustBoxSize(point_t &size, bool_point_t& frameRCfront, bool_point_t& frameRCend);

  /*!
  * The routine for creating wall \a Particle s
  */
  virtual void createParticles();
    
  /*!
  * Write out the final state
  */
  virtual ostream& write(ostream &s, int shift);
};

#endif
