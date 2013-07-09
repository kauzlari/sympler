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



#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"
#include "f_specific.h"
#include "colour_pair.h"
#include "pc_tube.h"

/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorTube>
		particle_creator_tube("ParticleCreatorTube");

#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation*) M_PHASE->parent())
#define PI M_PI

//---- Constructor/Destructor ----

ParticleCreatorTube::ParticleCreatorTube(Boundary *boundary) :
	ParticleCreatorFreePCalc(boundary) {
	init();
}

//---- Methods ----

void ParticleCreatorTube::adjustBoxSize(point_t &size,
		bool_point_t& frameRCfront, bool_point_t& frameRCend) {
	bool tube_defined = false, dd_defined = false;

	MSG_DEBUG("ParticleCreatorTube::adjustBoxSize", "size = " << size);

	if (m_radius != -1 && m_distance_rings != -1) {
		if (m_nring_particles < 1 || m_nrings < 1)
			throw gError("ParticleCreatorTube::adjustBoxSize", "At least one point is needed.");
		MSG_DEBUG("ParticleCreatorTube::adjustBoxSize", "m_ntube_points = (" << m_nring_particles << ", " << m_nrings << ")");

		tube_defined = true;
	}
	dd_defined = true;
}

void ParticleCreatorTube::setup() {
  ParticleCreatorFreePCalc::setup();

  double check_angle;
  // the rotation of the rings has to be a whole devidor of 180
  if (m_rotation_angle == 0) {
    m_rotation_fraction = 1;
  } else {
    check_angle = 360/m_rotation_angle;
    m_rotation_fraction = (int)(360/m_rotation_angle);
    if (check_angle != m_rotation_fraction)
      throw gError("ParticleCreatorTube::createParticles", "the rotation of the rings has to be a whole devidor of 180.");
  }
  m_rotation_angle = m_rotation_angle*PI/180;
  
  if (m_nrings < 3)
    throw gError("ParticleCreatorTube::createParticles", "Please put at least 3 rings.");
  if (m_radius <= 0)
    throw gError("ParticleCreatorTube::createParticles", "The tube radius must be larger than 0.");
  if (m_distance_rings <= 0)
    throw gError("ParticleCreatorTube::createParticles", "The distance between rings must be larger than 0.");
  
  if (m_connection_name == "")
    throw gError("ParticleCreatorTube::createParticles", "Please define connectorName.");
  if (m_triplet_force_name == "")
    throw gError("ParticleCreatorTube::createParticles", "Please define Angular force Name.");
  if (!((m_tube_configuration == "free") || (m_tube_configuration == "single clamped") || (m_tube_configuration == "double clamped") || (m_tube_configuration == "frozen")))
    throw gError("ParticleCreatorTube::createParticles", "Please define tube configuration: free, single clamped, double clamped or frozen.");
  
  ColourPair* cp = M_PHASE->manager()->cp(m_colour, m_colour);
  
  /*size_t listIndex =*/ cp -> createConnectedListIndexAndWrite(m_connection_name);
  /*size_t listIndex =*/ M_PHASE -> createTripletListIndexAndWrite(m_triplet_force_name);

}

void ParticleCreatorTube::createParticles() {
	ManagerCell *manager= M_PHASE->manager();
	RandomNumberGenerator rng;
	point_t offset, box_size;
	double angle;

	// randomize the velocities according to the input file
	if (!m_randomize) {
		rng.setSeed(m_seed);
	} else {
		rng.setSeed(time(0));
	}

	/* We need the manager to check which group a particle has. */

	initTransform();
	// creation of FREE particles
	MSG_DEBUG("ParticleCreatorTube::createParticles", "Creating tube particles.");

	//creating the first ring
	angle = PI*(m_nring_particles-2)/m_nring_particles;
	matrix_t vect(m_nring_particles, 3);
	matrix_t temp_vect(1, 2);
	// create the rotation angle
	matrix_t rotation_matrix(3, 3);
	createRotationMatrix(rotation_matrix);
	// first particle
	// it is ofset because we can't have negative positions here
	// the rotation vector for creatng the next particles
	vect(0, 0) = m_radius;
	vect(0, 1) = 0;
	vect(0, 2) = m_distance_rings/2;
	for (int j = 1; j < m_nring_particles; j++) {
		vect(j, 0) = vect(j-1, 0)*cos(PI-angle)-vect(j-1, 1)*sin(PI-angle);
		vect(j, 1) = vect(j-1, 1)*cos(PI-angle)+vect(j-1, 0)*sin(PI-angle);
		vect(j, 2) = m_distance_rings/2;
	}

	// creating the tube in the center of the box	
	offset = M_BOUNDARY->boundingBox().corner1;
	box_size = M_BOUNDARY->boundingBox().size();

	for (int i = 0; i < m_nrings; i++) {
	  for (int j = 0; j < m_nring_particles; j++) {
	    point_t pos = { { { vect(j, 0) + m_radius, vect(j, 1) + m_radius,
				vect(j, 2) + i*m_distance_rings } } };
	    
	    point_t pos_r = { { { rotation_matrix(0, 0)*pos.x+rotation_matrix(
									      0, 1)*pos.y+rotation_matrix(0, 2)*pos.z+box_size[0]/2,
				  rotation_matrix(1, 0)*pos.x+rotation_matrix(1, 1)*pos.y
				  +rotation_matrix(1, 2)*pos.z+box_size[1]/2,
				  rotation_matrix(2, 0)*pos.x+rotation_matrix(2, 1)*pos.y
				  +rotation_matrix(2, 2)*pos.z+box_size[2]/2 } } };
	    
	    Cell *c;
	    Particle p;
	    p.r = pos_r;
	    p.setColour(m_colour);
	    transformPos(p);
	    c = manager->findCell(p.r);
	    if (c) {
	      if (M_BOUNDARY->isInside(p.r)) {
		p.g = c->group();
		
		for (int w = 0; w < SPACE_DIMS; w++)
		  p.v[w] = rng.uniform() - 0.5;
		// Register particles
		if (m_tube_configuration =="free") {
		  m_particles[p.g].newEntry() = p;
		  if (i == 0)
		    MSG_DEBUG("ParticleCreatorTube::createParticles", "Creating free tube.");
		} else if (m_tube_configuration == "single clamped") {
		  if (i == 0) // the first ring is frozen
		    {
		      m_particles_frozen[p.g].newEntry() = p;
		      MSG_DEBUG("ParticleCreatorTube::createParticles", "Creating single clamped tube.");
		    }
		  if (i>0)
		    m_particles[p.g].newEntry() = p;
		  
		} else if (m_tube_configuration == "double clamped") {
		  if (i == 0) // the first ring is frozen
		    {
		      m_particles_frozen[p.g].newEntry() = p;
		      MSG_DEBUG("ParticleCreatorTube::createParticles", "Creating double clamped tube.");
		    }
		  if (i>0)
		    if (i == m_nrings-1) // the last ring is frozen too
		      m_particles_frozen[p.g].newEntry() = p;
		    else
		      m_particles[p.g].newEntry() = p;
		  
		} else if (m_tube_configuration == "frozen") {
		  m_particles_frozen[p.g].newEntry() = p;
		  MSG_DEBUG("ParticleCreatorTube::createParticles", "Creating frozen tube.");
		} else {
		  throw gError("ParticleCreatorTube::createParticles: This line should not be reaced.");
		}
	      } else {
		MSG_DEBUG("ParticleCreatorTube::", "M_BOUNDARY->isInside(p.r) not!"<<" "<< p.r[0]<<" "<<p.r[1]<<" "<<p.r[2]);
	      }
	    } else {
	      throw gError("ParticleCreatorTube::createParticles", "Particle not in boundaries!!!");
	    }
	  }
	  
	  // rotate the whole thing
	  temp_vect(0, 0) = vect(0, 0);
	  temp_vect(0, 1) = vect(0, 1);
	  vect(0, 0) = temp_vect(0, 0)*cos(m_rotation_angle)-temp_vect(0, 1)
	    *sin(m_rotation_angle);
	  vect(0, 1) = temp_vect(0, 1)*cos(m_rotation_angle)+temp_vect(0, 0)
	    *sin(m_rotation_angle);
	  for (int j = 1; j < m_nring_particles; j++) {
	    vect(j, 0) = vect(j-1, 0)*cos(PI-angle)-vect(j-1, 1)*sin(PI-angle);
	    vect(j, 1) = vect(j-1, 1)*cos(PI-angle)+vect(j-1, 0)*sin(PI-angle);
	  }
	}

	scaleVels();

	Particle *first_p, *first_frozen_p;
	flushParticles(&first_frozen_p, &first_p);
	// FIXME: the last two arguments are members, so unnecessary, no?
	connectTubeParticles(first_frozen_p, first_p, m_nring_particles, m_nrings);
	registerPulledParticles(first_p, m_nring_particles, m_nrings);
}

void ParticleCreatorTube::init() {
	m_properties.setClassName("ParticleCreatorTube");

	m_properties.setDescription("Generates particles that form a tube."
		" The tube can be either free, single clamped, double clamped or frozen."
		" If the tube is double clamped, a force can be applied in the center."
		" If the tube is single clamped, a force can be applied on the free tip."
		" Note that the size of the simulation box can't be modified by this"
		" ParticleCreator.");
	DOUBLEPC (temperature, m_temperature, 0,
			"Initial temperature of the particles. Note that this "
			"can be overriden by setting vel*.")
	;

	BOOLPC(randomize, m_randomize,
			"Set to 'no' to get the same random numbers for every simulation run.")
	;
	INTPC(seed, m_seed, 0,
			"The seed for the velocities. Only valid if randomize = no.")
	;
	INTPC(numberOfParticlesInRing , m_nring_particles, 0,
			"The number of particles in a ring of the tube")
	;

	INTPC(numberOfRingsInTube, m_nrings, 0,
			"The number of rings in a tube")
	;

	DOUBLEPC (radius, m_radius, 0,
			"Spacing of the particles in a ring.")
	;

	DOUBLEPC (distOfRings, m_distance_rings, 0,
			"Spacing of the rings of the tube.")
	;

	DOUBLEPC (ringRotationAngle, m_rotation_angle, -1,
			"The rotation of each ring in comparison to the former one. in degrees")
	;

	DOUBLEPC (rotateTheta, m_rotate_teta, -1,
			"The angle of rotation theta (angle to the z axis).")
	;

	DOUBLEPC (rotatePhi, m_rotate_phi, -1,
			"the angle of rotation phi (rotation on the xy plane).")
	;

	STRINGPC (connectorName, m_connection_name,
			"the name of the connection connecting the particles.")
	;

	STRINGPC (tripletForceName, m_triplet_force_name,
			"the name of the connection connecting the particles.")
	;

	STRINGPC (tubeConfiguration, m_tube_configuration,
			"The configuration of the tube: free, single clamped, double clamped or frozen.")
	;
	STRINGPC (pullForceName, m_pull_name,
			"The name of the force applied on the center ring of the tube if it is double clamped or last ring if it is single clamped.")
	;

	m_tube_configuration = "free";
	m_randomize = false;
	m_seed = 1;
	m_temperature = 1;
	m_rotate_phi = 0;
	m_rotate_teta = 0;
	m_nring_particles = 0;
	m_nrings = 0;
	m_radius = 0;
	m_distance_rings = 0;
	m_rotation_angle = 0;
}

void ParticleCreatorTube::scaleVels() {
	for (map<int, ParticleList>::iterator g = m_particles.begin(); g
			!= m_particles.end(); g++) {
		g->second.scaleVels(m_temperature);
	}
}

/* connectTubeparticles creates connections between particle pairs which match the   connections of a tube.
 it adds to the pairlist of the specified connection.
 Be careful: you can only connect after the particles are on the list (flushed)
 input: 
 the first particle in the tube
 the number of particles in a ring
 the number of rings
 */

void ParticleCreatorTube::connectTubeParticles(Particle *first_frozen_p,
		Particle *first_p, int nparticles_on_ring, int ntube_rings) {
  ColourPair* cp = M_PHASE->manager()->cp(m_colour, m_colour);

  if (m_tube_configuration == "free") {
    connectInRing(/*OLD*//*connector,*/cp, first_p, nparticles_on_ring, ntube_rings);
    connectNet(/*OLD*//*connector,*/cp, first_p, nparticles_on_ring, ntube_rings);
    addAngularForceFree(first_p, nparticles_on_ring, ntube_rings,
			m_rotation_fraction /*OLD*//*, tripletForce*/);
  }
  if (m_tube_configuration == "single clamped") {
    connectFrozenRing1(/*OLD*//*connector,*/cp, first_frozen_p, first_p,
		       nparticles_on_ring, ntube_rings);
    connectInRing(/*OLD*//*connector,*/cp, first_p, nparticles_on_ring, ntube_rings-1);
    connectNet(/*OLD*//*connector,*/cp, first_p, nparticles_on_ring, ntube_rings-1);
    addAngularForce1(first_frozen_p, first_p, nparticles_on_ring,
		     ntube_rings, m_rotation_fraction /*OLD*//*, tripletForce*/);
  }
  if (m_tube_configuration == "double clamped") {
    connectFrozenRing2(/*OLD*//*connector,*/cp, first_frozen_p, first_p,
		       nparticles_on_ring, ntube_rings);
    connectInRing(/*OLD*//*connector,*/cp, first_p, nparticles_on_ring, ntube_rings-2);
    connectNet(/*OLD*//*connector,*/cp, first_p, nparticles_on_ring, ntube_rings-2);
    addAngularForce2(first_frozen_p, first_p, nparticles_on_ring,
		     ntube_rings, m_rotation_fraction /*OLD*//*, tripletForce*/);
  }
  if (m_tube_configuration == "frozen") {
  }
}

void ParticleCreatorTube::connectNet(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *first_p, int nparticles_on_ring, int ntube_rings) {
  Particle *first_p_in_ring = first_p;
  Particle *first_ring_p = first_p_in_ring;
  Particle *second_ring_p = first_p_in_ring;
  
  // find the right list in the ColourPair to which to add the connections
  size_t listIndex = cp -> connectedListIndex(m_connection_name);
  
  for (int k = 0; k < nparticles_on_ring; k++)
    second_ring_p = second_ring_p->next;
  
  for (int i = 0; i < ntube_rings-1; i++) {
    // OLD: 1 line
    // 	  connector->addPairToConnectionAndWrite(first_ring_p, second_ring_p);
    // NEW: 1 line
    cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
//     MSG_DEBUG("ParticleCreatorTube::connectNet", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
    for (int k = 0; k < nparticles_on_ring-1; k++) {
      first_ring_p = first_ring_p->next;
      cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
//       MSG_DEBUG("ParticleCreatorTube::connectNet", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
      second_ring_p = second_ring_p->next;
      cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
//       MSG_DEBUG("ParticleCreatorTube::connectNet", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
    }
    cp->addPairToConnectionAndWrite(first_p_in_ring, second_ring_p, listIndex);
//     MSG_DEBUG("ParticleCreatorTube::connectNet", "Connected: " << first_p_in_ring->mySlot << "and "<< second_ring_p->mySlot);
    // finished connecting ring
    first_ring_p = first_ring_p->next;
    first_p_in_ring = first_ring_p;
    second_ring_p = second_ring_p->next;
  }
}

void ParticleCreatorTube::connectFrozenRing1(/*OLD*//*GenConnector *connector,*/ColourPair* cp,
		Particle *first_frozen_p, Particle *first_p, int nparticles_on_ring,
		int ntube_rings) {
	Particle *first_ring_p = first_frozen_p;
	Particle *second_ring_p = first_p;
	Particle *first_ring_p_static = first_frozen_p;

        // find the right list in the ColourPair to which to add the connections
        size_t listIndex = cp -> connectedListIndex(m_connection_name);

	/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
	for (int k = 0; k < nparticles_on_ring-1; k++) {
		first_ring_p = first_ring_p->next;
		/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 		MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
		second_ring_p = second_ring_p->next;
		/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 		MSG_DEBUG("ParticleCreatorTubE::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
	}
	/*connector*/cp->addPairToConnectionAndWrite(first_ring_p_static, second_ring_p, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p_static->mySlot << "and "<< second_ring_p->mySlot);
	// finished connecting first frozen ring
}

void ParticleCreatorTube::connectFrozenRing2(/*OLD*//*GenConnector *connector,*/ColourPair* cp,
		Particle *first_frozen_p, Particle *first_p, int nparticles_on_ring,
		int ntube_rings) {
	Particle *first_ring_p = first_frozen_p;
	Particle *second_ring_p = first_p;
	Particle *first_ring_p_static = first_frozen_p;

        // find the right list in the ColourPair to which to add the connections
        size_t listIndex = cp -> connectedListIndex(m_connection_name);

	/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
	for (int k = 0; k < nparticles_on_ring-1; k++) {
		first_ring_p = first_ring_p->next;
		/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 		MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
		second_ring_p = second_ring_p->next;
		/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 		MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
	}
	/*connector*/cp->addPairToConnectionAndWrite(first_ring_p_static, second_ring_p, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p_static->mySlot << "and "<< second_ring_p->mySlot);
	// finished connecting first frozen ring

	second_ring_p = first_frozen_p;
	for (int k = 0; k < nparticles_on_ring; k++) {
		second_ring_p = second_ring_p->next;
	}
	first_ring_p = first_p;
	for (int k = 0; k < (ntube_rings-3)*nparticles_on_ring; k++) {
		first_ring_p = first_ring_p->next;
	}
	first_ring_p_static = first_ring_p;
	/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
	for (int k = 0; k < nparticles_on_ring-1; k++) {
		first_ring_p = first_ring_p->next;
		/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
// 		MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
		second_ring_p = second_ring_p->next;
		/*connector*/cp->addPairToConnectionAndWrite(first_ring_p, second_ring_p, listIndex);
		MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p->mySlot << "and "<< second_ring_p->mySlot);
	}
	/*connector*/cp->addPairToConnectionAndWrite(first_ring_p_static, second_ring_p, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectNetFrozen", "Connected: " << first_ring_p_static->mySlot << "and "<< second_ring_p->mySlot);
	// finished connecting last frozen ring
}

void ParticleCreatorTube::connectInRing(/*OLD*//*GenConnector *connector,*/ColourPair* cp,
		Particle *first_p, int nparticles_on_ring, int ntube_rings) {
	Particle *first_p_in_ring = first_p;

	for (int i = 0; i < ntube_rings; i++) {
	  connectCyclic(/*connector*/cp, first_p_in_ring, nparticles_on_ring);
		for (int j = 0; j < nparticles_on_ring; j++)
			first_p_in_ring = first_p_in_ring->next;
	}
}

void ParticleCreatorTube::connectCyclic(/*OLD*//*GenConnector *connector,*/ColourPair* cp, Particle *p,
		int num) {
	Particle *first = p;

        // find the right list in the ColourPair to which to add the connections
        size_t listIndex = cp -> connectedListIndex(m_connection_name);

	for (int k = 0; k < num-1; k++) {
	  /*connector*/cp->addPairToConnectionAndWrite(p, p->next, listIndex);
// 		MSG_DEBUG("ParticleCreatorTube::connectCyclic", "Connected: " << p->mySlot << "and "<< p->next->mySlot);
		p = p->next;
	}
	/*connector*/cp->addPairToConnectionAndWrite(p, first, listIndex);
// 	MSG_DEBUG("ParticleCreatorTube::connectCyclic", "Connected: " << p->mySlot << "and "<< first->mySlot);
}

void ParticleCreatorTube::addAngularForceFree(Particle *first_p,
		int nparticles_on_ring, int ntube_rings, int rotation_fraction
					      /*OLD*/ /*, GenTriplet* triplet*/) {

        // find the right list in the Phase to which to add the triplets
        size_t listIndex = M_PHASE -> tripletListIndex(m_triplet_force_name);

	// this is a little complicated algorithm...
	int full_circle = rotation_fraction;
	Particle *first_ring_p, *p1;
	Particle *second_ring_p, *p2;
	Particle *third_ring_p, *p3;
	first_ring_p = first_p;
	second_ring_p = first_p;
	third_ring_p = first_p;

	for (int i = 0; i < full_circle; i++) {
		second_ring_p = second_ring_p->next;
	}
	third_ring_p = second_ring_p;
	for (int i = 0; i < full_circle; i++) {
		third_ring_p = third_ring_p->next;
	}

	for (int j = 0; j < ntube_rings-2*rotation_fraction/nparticles_on_ring; j++) {
		p1 = first_ring_p->next->next;
		p2 = second_ring_p->next;
		p3 = third_ring_p;
		/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
		for (int i = 0; i < nparticles_on_ring-3; i++) {
			p1 = p1->next;
			p2 = p2->next;
			p3 = p3->next;
			/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
		}
		// add the triplets of the first two in the ring
		p2 = p2->next;
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p, p2, p3, listIndex);
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p->next, second_ring_p, p3, listIndex);

		// advance the pointers only if this is not the last triplet
		if (j < ntube_rings-2*rotation_fraction/nparticles_on_ring-1) {
			first_ring_p = p1->next;
			second_ring_p = p2->next;
			third_ring_p = p3->next;
		}
	}
}

void ParticleCreatorTube::addAngularForce2(Particle *first_frozen_p,
		Particle *first_p, int nparticles_on_ring, int ntube_rings,
					   int rotation_fraction/*, GenTriplet* triplet*/) {

        // find the right list in the Phase to which to add the triplets
        size_t listIndex = M_PHASE -> tripletListIndex(m_triplet_force_name);

	// this is a little complicated algorithm...
	int full_circle = rotation_fraction;
	Particle *first_ring_p, *p1;
	Particle *second_ring_p, *p2;
	Particle *third_ring_p, *p3;
	first_ring_p = first_p;
	second_ring_p = first_p;
	third_ring_p = first_p;

	// FIXME: !!! this needs to be fixed for short tubes !!!

	// connect the frozen ring to the free particles	
	for (int i = 0; i < full_circle - nparticles_on_ring; i++) {
		second_ring_p = second_ring_p->next;
	}
	third_ring_p = second_ring_p;
	for (int i = 0; i < full_circle; i++) {
		third_ring_p = third_ring_p->next;
	}

	p1 = first_frozen_p->next->next;
	p2 = second_ring_p->next;
	p3 = third_ring_p;
	/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
	for (int i = 0; i < nparticles_on_ring-3; i++) {
		p1 = p1->next;
		p2 = p2->next;
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
	}
	// add the triplets of the first two in the ring
	p2 = p2->next;
	p3 = p3->next;
	/*triplet*/M_PHASE->addTripletAndWrite(first_frozen_p, p2, p3, listIndex);
	p3 = p3->next;
	/*triplet*/M_PHASE->addTripletAndWrite(first_frozen_p->next, second_ring_p, p3, listIndex);

	// connect the rest of the rings
	first_ring_p = first_p;
	second_ring_p = p2->next;
	third_ring_p = p3->next;

	for (int j = 0; j < ntube_rings-2-2*rotation_fraction/nparticles_on_ring; j++) {
		p1 = first_ring_p->next->next;
		p2 = second_ring_p->next;
		p3 = third_ring_p;
		/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
		for (int i = 0; i < nparticles_on_ring-3; i++) {
			p1 = p1->next;
			p2 = p2->next;
			p3 = p3->next;
			/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
		}
		// add the triplets of the first two in the ring
		p2 = p2->next;
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p, p2, p3, listIndex);
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p->next, second_ring_p, p3, listIndex);

		first_ring_p = p1->next;
		second_ring_p = p2->next;
		third_ring_p = p3->next;
	}
	// connect to the last ring which is frozen
	third_ring_p = first_frozen_p;
	for (int k = 0; k < nparticles_on_ring; k++) {
		third_ring_p = third_ring_p->next;
	}
	p1 = first_ring_p->next->next;
	p2 = second_ring_p->next;
	p3 = third_ring_p;
	/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
	for (int i = 0; i < nparticles_on_ring-3; i++) {
		p1 = p1->next;
		p2 = p2->next;
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
	}
	// add the triplets of the first two in the ring
	p2 = p2->next;
	p3 = p3->next;
	/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p, p2, p3, listIndex);
	p3 = p3->next;
	/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p->next, second_ring_p, p3, listIndex);

}
void ParticleCreatorTube::addAngularForce1(Particle *first_frozen_p,
		Particle *first_p, int nparticles_on_ring, int ntube_rings,
					   int rotation_fraction/*, GenTriplet* triplet*/) {

        // find the right list in the Phase to which to add the triplets
        size_t listIndex = M_PHASE -> tripletListIndex(m_triplet_force_name);

	// this is a little complicated algorithm...
	int full_circle = rotation_fraction;
	Particle *first_ring_p, *p1;
	Particle *second_ring_p, *p2;
	Particle *third_ring_p, *p3;
	first_ring_p = first_p;
	second_ring_p = first_p;
	third_ring_p = first_p;

	// connect the frozen ring to the free particles	
	for (int i = 0; i < full_circle - nparticles_on_ring; i++) {
		second_ring_p = second_ring_p->next;
	}
	third_ring_p = second_ring_p;
	for (int i = 0; i < full_circle; i++) {
		third_ring_p = third_ring_p->next;
	}

	p1 = first_frozen_p->next->next;
	p2 = second_ring_p->next;
	p3 = third_ring_p;
	/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
	for (int i = 0; i < nparticles_on_ring-3; i++) {
		p1 = p1->next;
		p2 = p2->next;
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
	}
	// add the triplets of the first two in the ring
	p2 = p2->next;
	p3 = p3->next;
	/*triplet*/M_PHASE->addTripletAndWrite(first_frozen_p, p2, p3, listIndex);
	p3 = p3->next;
	/*triplet*/M_PHASE->addTripletAndWrite(first_frozen_p->next, second_ring_p, p3, listIndex);

	// connect the rest of the rings
	first_ring_p = first_p;
	second_ring_p = p2->next;
	third_ring_p = p3->next;

	for (int j = 0; j < ntube_rings-1-2*rotation_fraction/nparticles_on_ring; j++) {
		p1 = first_ring_p->next->next;
		p2 = second_ring_p->next;
		p3 = third_ring_p;
		/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
		for (int i = 0; i < nparticles_on_ring-3; i++) {
			p1 = p1->next;
			p2 = p2->next;
			p3 = p3->next;
			/*triplet*/M_PHASE->addTripletAndWrite(p1, p2, p3, listIndex);
		}
		// add the triplets of the first two in the ring
		p2 = p2->next;
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p, p2, p3, listIndex);
		p3 = p3->next;
		/*triplet*/M_PHASE->addTripletAndWrite(first_ring_p->next, second_ring_p, p3, listIndex);

		// advance the pointers only if this is not the last triplet
		if (j < ntube_rings-2*rotation_fraction/nparticles_on_ring-1) {
			first_ring_p = p1->next;
			second_ring_p = p2->next;
			third_ring_p = p3->next;
		}
	}

}

void ParticleCreatorTube::registerPulledParticles(Particle *first_p, int nparticles_on_ring, int nrings) 
{
	Particle *p = first_p;
	Fspecific* force;
	
	if (m_tube_configuration == "single clamped") {
		force = (Fspecific*)M_SIMULATION->findForceFunction(m_pull_name);
		// get the first particle of the last ring
		for (int i = 0; i < nrings-2; i++)
			for (int j = 0; j < nparticles_on_ring; j++)
				p = p->next;
		for (int j = 0; j < nparticles_on_ring; j++) {
			force->addParticleToForce(p);
			MSG_DEBUG("ParticleCreatorTube::registerPull", "registered: " << p->mySlot);
			p = p->next;
		}
	}
	if (m_tube_configuration == "double clamped") {
		force = (Fspecific*)M_SIMULATION->findForceFunction(m_pull_name);
		// get the first particle of the center ring
		for (int i = 0; i < round(nrings/2)-1; i++)
			for (int j = 0; j < nparticles_on_ring; j++)
				p = p->next;
		for (int j = 0; j < nparticles_on_ring; j++) {
			force->addParticleToForce(p);
			MSG_DEBUG("ParticleCreatorTube::registerPull", "registered: " << p->mySlot);
			p = p->next;
		}
	}
}

void ParticleCreatorTube::flushParticles(Particle** first_frozen_p,
		Particle** first_p) {
	Phase *phase= M_PHASE;
	size_t counter = 0;
	if (!m_particles_frozen.empty()) {
	  for (map<int, ParticleList>::iterator g = m_particles_frozen.begin(); g
		 != m_particles_frozen.end(); g++) {
	    SL_FOR_EACH
	      (Particle, g->second,
	       transformVel(*__iSLFE);
	       if (__iSLFE == m_particles_frozen.begin()->second.first())
		 {
		   *first_frozen_p = phase->addFrozenParticle(*__iSLFE);
		 }
	       else
		 phase->addFrozenParticle(*__iSLFE);
	       ++counter;
	       )
	      ;
	  }
	}

	if (!m_particles.empty()) {
	  for (map<int, ParticleList>::iterator g = m_particles.begin(); g
		 != m_particles.end(); g++) {
	    SL_FOR_EACH
	      (Particle, g->second,
	       transformVel(*__iSLFE);
	       if (__iSLFE == m_particles.begin()->second.first())
		 {
		   *first_p = phase->addParticle(*__iSLFE);
		 }
	       else
		 phase->addParticle(*__iSLFE);
	       ++counter;
	       )
	      ;
	  }
	}

	MSG_DEBUG("ParticleCreatorTube::flushParticles for " << m_properties.name(), counter << " particles added");
	m_particles_frozen.clear();
	m_particles.clear();
	//MSG_DEBUG("ParticleCreatorTube::flushParticles", "after m_particles.clear()");
}

void ParticleCreatorTube::createRotationMatrix(matrix_t& matrix) {
	m_rotate_phi = m_rotate_phi/180*PI;
	m_rotate_teta = m_rotate_teta/180*PI;
	matrix(0, 0) = cos(m_rotate_phi)*cos(m_rotate_teta);
	matrix(0, 1) = sin(m_rotate_phi);
	matrix(0, 2) = -cos(m_rotate_phi)*sin(m_rotate_teta);;
	matrix(1, 0) = -cos(m_rotate_teta)*sin(m_rotate_phi);
	matrix(1, 1) = cos(m_rotate_phi);
	matrix(1, 2) = sin(m_rotate_teta)*sin(m_rotate_phi);
	matrix(2, 0) = sin(m_rotate_teta);
	matrix(2, 1) = 0;
	matrix(2, 2) = cos(m_rotate_teta);
}
