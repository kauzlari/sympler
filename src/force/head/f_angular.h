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



#ifndef _F_ANGULAR_H
#define _F_ANGULAR_H

#include "gen_triplet.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


/*!
 * A force between three particles. depends on the angle of the triplet F=k(theta-thetaEq)
 */
class FAngular : public GenTriplet
{
	protected:
/*!
 * Initialize the property list
 */
	void init();
/*!
 * m_k angular stiffnes
*/
	double m_k;
/*!
 * m_thetaEq equilibrium angle of the angular spring
*/
	double m_thetaEq;
/*!
 * m_cosEq Cosine of the equilibrium angle 
*/
	double m_cosEq;

/*!
 * m_periodic use periodic boundarys. (03.08.2015) Cannot distiguish directions yet (implementation)
*/	
	bool m_periodic;

	public:
/*!
 * Constructor
 */
	FAngular(Simulation *simulation);
	FAngular();

/*!
 * Destructor
 */
	virtual ~FAngular();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no) {}
#endif

/*!
 * Computes the force on the three particles.
 * @param force_index The index within the force history
 */
  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

/*!
 * Setup shortly before simulation starts
 */
	virtual void setup();

/*!
 * Setup assigning the right triplet list to this force
 */
/* virtual void setupAfterParticleCreation(); */


};

#endif
