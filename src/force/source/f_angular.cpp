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




#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "triplet.h"
#include "f_angular.h"

// #include "valgrind/memcheck.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

#define PI M_PI

const GenFTypeConcr<FAngular> angular_force("FAngular");

//---- Constructors/Destructor ----

FAngular::FAngular()
{
	throw gError("FAngular::FAngular(default)", "Should not be called. Contact the programmer.");
}

FAngular::FAngular(Simulation *simulation): GenTriplet(simulation)
{
	init();
}


FAngular::~FAngular()
{
}


void FAngular::init()
{
	m_properties.setClassName("FAngular");

	m_properties.setDescription(
			"This is an angular force: F=k(cos(theta)-cos(theta_eq)).");

	DOUBLEPC(K, m_k, -1,"Spring force constant.");

	DOUBLEPC(thetaEq, m_thetaEq, -1,"equilibrium angle (in degrees).");

	BOOLPC(periodic, m_periodic, "Should periodic boundary conditions be applied to the connections?");

	/*!
	 * m_k angular stiffnes
	*/
	m_k = 0;

	/*!
	 * m_thetaEq equilibrium angle of the angular spring
	*/	
	m_thetaEq = 0;

	/*!
	 * m_periodic use periodic boundarys. (03.08.2015) Cannot distiguish directions yet (implementation)
	*/
	m_periodic = true; 

	/*!
	 * This Force is neither a particle force nor a pair force!
	*/
	m_is_pair_force = false;
	m_is_particle_force = false;
}


//---- Methods ----

void FAngular::computeForces(int force_index)
{
  //   MSG_DEBUG("FAngular::computeForces", "START");
  // FIXME: not so nice to do that here, but currently the only possibility. See also GenTriplet::setupAfterParticleCreation(). Is the problem solved if this function will be called for each triplet (when parallelised)?  
  if(!m_TripletList) {
    MSG_DEBUG("FAngular::computeForces", "requesting triplet list");
    m_TripletList = M_PHASE -> returnTripletList(m_force_name);
    MSG_DEBUG("FAngular::computeForces", "got triplet list");
}

  point_t boxSize = M_PHASE->boundary()->boundingBox().size();
//   double size;

// Old: 1 line
//   for (tripletListItr p = m_TripletList.begin(); p != m_TripletList.end(); ++p) 
  for (tripletListItr p = m_TripletList->begin(); p != m_TripletList->end(); ++p) 
    {	
      point_t b1, b2; // vector
      double c11 = 0;
      double c12 = 0;
      double c22 = 0; // scalar product
//       double abs_b1, abs_b2;
      double cos_a; //the cosine of the angle
      double F;
//       double a;
      point_t force_a, force_b, force_c;
      for (int _i = 0; _i < SPACE_DIMS; _i++)
	{
	  b1[_i] = p->b->r[_i] - p->a-> r[_i];
	  b2[_i] = p->c->r[_i] - p->b-> r[_i];
	  // periodic BCs
	  if(m_periodic == true)
	    {
	      if(b1[_i] > 0.5*boxSize[_i]) b1[_i] -= boxSize[_i]; 
	      if(b1[_i] < -0.5*boxSize[_i]) b1[_i] += boxSize[_i]; 
	      if(b2[_i] > 0.5*boxSize[_i]) b2[_i] -= boxSize[_i]; 
	      if(b2[_i] < -0.5*boxSize[_i]) b2[_i] += boxSize[_i]; 
	  }

	  c11 += b1[_i] * b1[_i];
	  c12 += b1[_i] * b2[_i];
	  c22 += b2[_i] * b2[_i];
	}
      double invAbsC11c22 = 1/sqrt(c11*c22);

      cos_a = c12*invAbsC11c22;

//        MSG_DEBUG("FAngular::calculate force", "particle " << p->b->mySlot << " cos(theta) = " << cos_a);
//       MSG_DEBUG("FAngular::calculate force", "particle " << p->b->mySlot << " m_thetaEq = " << m_thetaEq);
//       MSG_DEBUG("FAngular::calculate force", "particle " << p->b->mySlot << " cosEQ = " << cos(PI*m_thetaEq/180));
//        MSG_DEBUG("FAngular::calculate force", "particle " << p->b->mySlot << " m_cosEq = " << m_cosEq);
      
      
      F = -m_k*(cos_a-m_cosEq);
      
      
      force_a = F*(c12/c11*b1-b2)*invAbsC11c22;
//        MSG_DEBUG("FAngular::calculate force", "particle " << p->a->mySlot << " force_a " << force_a);
      force_c = F*(b1-c12/c22*b2)*invAbsC11c22;
//        MSG_DEBUG("FAngular::calculate force", "particle " << p->c->mySlot << " force_c " << force_c);
      force_b = -1*force_a - force_c;
//        MSG_DEBUG("FAngular::calculate force", "particle " << p->b->mySlot << " force_b " << force_b);
      p->a->force[force_index] += force_a;
      p->b->force[force_index] += force_b;
      p->c->force[force_index] += force_c;
    }
}


#ifndef _OPENMP
void FAngular::computeForces(Pairdist* pair, int force_index)
#else
void FAngular::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FAngular::computeForces", "Fatal error: do not call FAngular::computeForces(Pairdist* pair, int force_index)!!! Please contact the programmer!");
}


void FAngular::computeForces(Particle* part, int force_index)
{
  throw gError("FAngular::computeForces", "Fatal error: do not call FAngular::computeForces(Particle* part, int force_index)!!! Please contact the programmer!");
}


void FAngular::setup()
{
  GenTriplet::setup();
  m_cosEq = cos(PI-PI*m_thetaEq/180.);	
  // 	MSG_DEBUG("FAngular::setup", "m_cosEqDiff" << m_cosEq-cos(PI*m_thetaEq/180.));
  // 	MSG_DEBUG("FAngular::setup", "m_cosEqDiff2" << cos(PI*m_thetaEq/180)-cos(PI*m_thetaEq/180));
}

