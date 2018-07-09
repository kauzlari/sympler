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


#include "f_particle.h"
#include "f_specific.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<Fspecific> fspecific("Fspecific");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

Fspecific::Fspecific(Simulation *simulation) :
	GenF(simulation) {
	init();
}

Fspecific::~Fspecific() {
}

void Fspecific::computeForces(int force_index) {
	for (SpecificParticleListItr p = m_SpecificParticleList.begin(); p
			!= m_SpecificParticleList.end(); ++p)
		(*p)->force[force_index] += m_force;
}


void Fspecific::computeForces(Particle* part, int force_index)
{
  throw gError("Fspecific::computeForces", "Fatal error: do not call Fspecific::computeForces(Particle* part, int force_index)!!! Please contact the programmer!");
}


#ifndef _OPENMP
void Fspecific::computeForces(Pairdist* pair, int force_index)
#else
void Fspecific::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("Fspecific::computeForces", "Fatal error: do not call Fspecific::computeForces(Pairdist* pair, int force_index)!!! Please contact the programmer!");
}


void Fspecific::init()
{
	m_properties.setClassName("Fspecific");
	m_properties.setDescription("Constant one-particle force.works on specific particles and not on a species");

	POINTPC
	(forceField, m_force,
			"Constant force field.")
	;

	m_properties.setName("Fspecific");

	STRINGPC
	(species, m_species," ")
	;

	m_species = "UNDEF";

	m_force.assign(0);

  m_is_pair_force = false;
  m_is_particle_force = false;
}

void Fspecific::setup() {
	GenF::setup();
}

/* add a particle to the list for the force field
 */
void Fspecific::addParticleToForce(Particle *p) {
	ofstream file_stream;
	string file_name;
	file_name = M_SIMULATION->name();
	// OLD: 2 lines; what is it good for?
// 	file_name += "_";
// 	file_name += M_MANAGER->species(p->c);
	file_name += "_connector.con";
	file_stream.open(file_name.c_str(), ios_base::app);
	file_stream.seekp(0, std::ios_base::end);
	file_stream << "single " << m_force_name << " " << p -> mySlot << " ";
	if (p->isFrozen)
		file_stream << "frozen" << " " << endl;
	else
		file_stream << "free" << " " << endl;

	file_stream.close();

	m_SpecificParticleList.push_back(p);

	MSG_DEBUG("Connector::addParticle to force field","");
}
