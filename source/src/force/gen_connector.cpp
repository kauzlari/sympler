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
#include "pairdist.h"
#include "pair_list.h"
#include "controller.h"
#include "manager_cell.h"
// #include "integrator_position.h"
// #include "integrator_velocity_verlet.h"
// #include "integrator_static.h"
#include "gen_connector.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

/*---- Class Connector ----*/

GenConnector::GenConnector() {
	throw gError("GenConnector::GenConnector(default)", "Should not be called. Contact the programmer.");
}

GenConnector::GenConnector(Simulation *simulation) :
	GenF(simulation), m_cp(NULL), m_connectedList(NULL) {
	init();
}

/*
 GenConnector::GenConnector(Simulation *simulation, ColourPair* cp) : GenF(simulation), m_cp(cp), m_connectedList(new PairList(cp))
 {
 init();
 }
 */
GenConnector::~GenConnector() {
  //OLD: 1 line; not responsible anymore
// 	delete(m_connectedList);
}

void GenConnector::init() {
	m_properties.setName("genConnector");
	m_properties.setDescription("Base class for all connection forces.");

	STRINGPC
	(species1, m_species.first,
			"First species, this force should act on.")
	;

	STRINGPC
	(species2, m_species.second,
			"Second species, this force should act on.")
	;

	m_species.first = "UNDEF";
	m_species.second = "UNDEF";
}

void GenConnector::setup() {
	GenF::setup();

	if (m_species.first == "UNDEF")
		throw gError("GenConnector::setup", "'species1' undefined in " + m_properties.name());
	if (m_species.second == "UNDEF")
		throw gError("GenConnector::setup", "'species2' undefined in " + m_properties.name());

	m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

}

void GenConnector::setupAfterParticleCreation() {

  // FIXME: The following doesn't work (but would be nice) because this function is called BEFORE the corresponding one from, e.g., a PCConnector.

// MSG_DEBUG("GenConnector::setupAfterParticleCreation", "before: m_connectedList = " << m_connectedList);
  
// // the following works as long as the setup of this force is called after the one of the modules creating the connected lists
//   m_connectedList = m_cp->connectedList(m_force_name);

// MSG_DEBUG("GenConnector::setupAfterParticleCreation", "after: m_connectedList = " << m_connectedList);
}



