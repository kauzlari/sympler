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



using namespace std;

#include <algorithm>

#include "phase.h"
#include "controller.h"
#include "manager_cell.h"
#include "integrator_position.h"
#include "gen_triplet.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

/*---- Class triplet force ----*/

GenTriplet::GenTriplet() {
	throw gError("GenTriplet::Triplet(default)", "Should not be called. Contact the programmer.");
}

GenTriplet::GenTriplet(Simulation *simulation) :GenF(simulation), m_TripletList(NULL) 
{
	init();
}

GenTriplet::~GenTriplet() {
}

void GenTriplet::init() {
	m_properties.setName("genTriplet");
	m_properties.setDescription("Base class for all triplet forces.");

  m_is_pair_force = false;
  m_is_particle_force = false;
}

void GenTriplet::setup() {
	GenF::setup();
// 	m_TripletList = M_PHASE -> returnTripletList(m_force_name);

// reason for commenting this out see GenTriplet::setupAfterParticleCreation()    
//   M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

}

void GenTriplet::setupAfterParticleCreation() {

  // FIXME: The following doesn't work (but would be nice) because this function is called BEFORE the corresponding one from, e.g., a PCConnector.

//   m_TripletList = M_PHASE -> returnTripletList(m_force_name);
}
