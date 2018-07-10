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


	
#ifndef __PC_TUBE_H
#define __PC_TUBE_H

#include "phase.h"
#include "pc_free_pcalc.h"
#include "boundary.h"	
	

//---- Classes ----

class ParticleCreatorTube: public ParticleCreatorFreePCalc
{
 protected:
  /* list of frozen particles */
  map<int, ParticleList> m_particles_frozen;
  /* Number of particles in each ring */
  int m_nring_particles;
  /* Number of rings in the tube */
  int m_nrings;
  /* distance between points in a ring */
  double m_radius;
  /* distance between rings */
  double m_distance_rings;
  /* ring rotation fraction */
  double m_rotation_angle;
  int m_rotation_fraction;
  /* tube angle */
  double m_rotate_teta, m_rotate_phi;
  
  double m_temperature;
  
  /* If randomize = false then seed = 1 (for rand) */
  bool m_randomize;
  int m_seed;	
  
  /* If force_boundary_size is on, the boundary will be resized to
     fit the specified information. */
  bool m_force_boundary_size;
  
  string m_connection_name;
  string m_triplet_force_name;
  string m_pull_name;
  string m_tube_configuration;
  
  void init();
  
  virtual void scaleVels();
  
  void connectTubeParticles(Particle *first_frozen_p,Particle *first_p,
			    int nparticles_on_ring,
			    int ntube_rings);
  void connectNet(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *first_p,int nparticles_on_ring, int ntube_rings);
  void connectInRing(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *first_p,int nparticles_on_ring, int ntube_rings);
  
  void connectFrozenRing1(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *first_frozen_p, Particle *first_p,int nparticles_on_ring, int ntube_rings);
  void connectFrozenRing2(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *first_frozen_p, Particle *first_p,int nparticles_on_ring, int ntube_rings);
  void connectCyclic(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *p,int num);
  
  void addAngularForceFree(Particle *first_p,
			   int nparticles_on_ring, int ntube_rings,int m_rotation_fraction
			   /*OLD*//*,GenTriplet *tripletForce*/);
  void addAngularForce1(Particle *first_frozen_p, Particle *first_p,
			int nparticles_on_ring, int ntube_rings,int m_rotation_fraction
			/*OLD*//*,GenTriplet *tripletForce*/);
  void addAngularForce2(Particle *first_frozen_p, Particle *first_p,
			int nparticles_on_ring, int ntube_rings,int m_rotation_fraction
			/*OLD*//*,GenTriplet *tripletForce*/);
  
  void registerPulledParticles(Particle *first_p,int nparticles_on_ring, int ntube_rings);
  void createRotationMatrix(matrix_t& matrix);
  
 public:
  ParticleCreatorTube(Boundary *boundary);
  
  virtual void adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend);
  
  virtual void createParticles();
  virtual void flushParticles(Particle** first_frozen_p, Particle** first_p);
  
  virtual int nOfFreeP() {
    return m_nring_particles*m_nrings;
  }   
  
  virtual void setup();
};

#endif

