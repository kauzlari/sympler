/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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
#include "gen_quintet.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

/*---- Class quintet force ----*/

GenQuintet::GenQuintet() {
	throw gError("GenQuintet::Quintet(default)", "Should not be called. Contact the programmer.");
}

GenQuintet::GenQuintet(Simulation *simulation) :GenF(simulation), m_QuintetList(NULL) 
{
	init();
}

GenQuintet::~GenQuintet() {
}

void GenQuintet::init() {
	m_properties.setName("genQuintet");
	m_properties.setDescription("Base class for all quintet forces.");

  m_is_pair_force = false;
  m_is_particle_force = false;
}

void GenQuintet::setup() {
	GenF::setup();

// reason for commenting this out see GenQuintet::setupAfterParticleCreation()    
//   M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

}

void GenQuintet::setupAfterParticleCreation() {

  // FIXME: The following doesn't work (but would be nice) because this function is called BEFORE the corresponding one from, e.g., a PCConnector.

}
