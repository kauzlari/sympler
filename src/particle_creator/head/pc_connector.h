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




/// This file contains the ParticleConnectorFile class,
/// a particle creator which connects existing particles
/// by specific, spring and triplet forces.
/// The particles are previously defined by any of the
/// other particles creators (almost). ParticleConnectorFile
/// then assigned forces (according to a file) to those
/// particles. So strictly speaking, it does not create any
/// particle at all.

#ifndef PC_CONNECTOR_H_
#define PC_CONNECTOR_H_

#include <vector>
#include <iostream>
#include "boundary.h"
#include "particle_creator.h"
#include "particle_list.h"

//---- Classes ----

/////////////////////////////////////////////////////////////////////
///  \a ParticleConnectorFile reads information about particle
///  connections from a file and applies the corresponding forces to
///  particles created by another \a ParticleCreator. So strictly
///  speaking, it does not create any particle at all.
/////////////////////////////////////////////////////////////////////
class ParticleConnectorFile : public ParticleCreator {
protected:
	///
	/// Initialize the property list
	///
	void init();

	///
	/// Vector of \a ParticleList s for the given species
	///	
	vector<ParticleList*> m_particlesToConnect;
	///
	/// Vector of \a ParticleList s for the given species
	///	
	vector<ParticleList*> m_frozenParticlesToConnect;

	///
	/// Get a particle from \a particlesToConnect
	/// @param n The index of the particle (starting from 0
	///   and spanning all particle lists).
	///	
	Particle* getParticleFromNumber(size_t n, string freeOrFrozen);

public:

	///
	/// Constructor
	/// Cannot be used.
	///
	ParticleConnectorFile();

	///
	/// Constructor
	/// @param boundary The \a Boundary this \a ParticleCreator belongs to
	///
	ParticleConnectorFile(Boundary*);

	///
	/// Destructor
	///
	virtual ~ParticleConnectorFile();

	///
	/// Overwrite old setup to remove check for species
	///
	virtual void setup();

	///
	/// Create (no) initial particles.
	///
	virtual void createParticles() {
	}

	///
	/// Setup connector forces
	///
	virtual void setupAfterParticleCreation();

	///
	/// Put setup into restart file
	/// Empty here, since done by \a ParticleCreatorFree.
	///
	virtual ostream &write(ostream &s, int shift = 0) {
		return s;
	}

	///
	/// File with connection data.
	///
	string m_filename;

};

#endif /*PC_CONNECTOR_H_*/
