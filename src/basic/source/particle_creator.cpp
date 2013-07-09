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



#include <fstream>
#include <iostream>

#include "phase.h"
#include "manager_cell.h"
// #include "simulation.h"

#include "particle_creator.h"


#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())
#define M_MANAGER M_PHASE->manager()

bool ParticleCreator::s_createFrozenParts;

REGISTER_SMART_ENUM
(ParticleCreator_Factory,
 "ParticleCreators are in charge of placing particles into the simulation box before the simulation "
 "starts. Most ParticleCreators create so called \"free\" particles, in the sense that all the particles' degrees of freedom are allowed to change. Some ParticleCreators can also create \"frozen\" particles. For frozen particles, all particle properties are \"frozen\" to the initial values and will never change. Frozen particles may be useful for setting certain kinds of boundary conditions."
);


ParticleCreator::ParticleCreator()
    : Node()
{
    MSG_DEBUG("ParticleCreator::ParticleCreator", "Shouldn't be called! FATAL!");
    // This should never be called!
	abort();
}

ParticleCreator::ParticleCreator(Boundary *boundary)
    : Node((Node*) boundary)
{
	init();
}



ParticleCreator::~ParticleCreator()
{
}



//--- Methods ---



void ParticleCreator::init()
{
	m_properties.setClassName("ParticleCreator");
	STRINGPC(species, m_species, "Defines the species, the ParticleCreator should create."
           " See the general description of this ParticleCreator for possible special"
           " functionalities concerning this attribute.");
	m_species = "UNDEF";
}


void ParticleCreator::setup()
{
	if(m_species == "UNDEF")
		throw gError("You must define 'species' for " + m_properties.name());

  m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
  s_createFrozenParts = false;
}


void ParticleCreator::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend)
{
}


void ParticleCreator::flushParticles()
{
#if 0
    Phase *phase = M_PHASE;

    for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
    //    for (map<int, list<Particle> >::iterator g = m_particles.begin(); g != m_particles.end(); g++){ 
        SL_FOR_EACH
            (Particle,
             g->second,
             if (M_BOUNDARY->isInside(__iSLFE->r))
                 phase->addParticle(*__iSLFE);
             );
    }

    /* Make sure we don't use too much memory while the simulation
       is running. */
    m_particles.clear();
#endif
}


void ParticleCreator::createMoreParticlesAtStart(size_t n)
{
  for (size_t i = 0; i < n; i++)
    createParticle(false);
}


