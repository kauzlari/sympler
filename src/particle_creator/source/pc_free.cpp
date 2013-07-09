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


#include "pc_free.h"

#include "simulation.h"
#include "manager_cell.h"

bool ParticleCreatorFree::toFileDone = false;


#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation*) M_PHASE->parent())
#define M_MANAGER M_PHASE->manager()

ParticleCreatorFree::ParticleCreatorFree()
	: ParticleCreator()
{
  throw gError("ParticleCreatorFree::ParticleCreatorFree()", "Forbidden call. Contact programmers.");
}

ParticleCreatorFree::ParticleCreatorFree(Boundary *boundary)
	: ParticleCreator(boundary)
{
	init();
}


ParticleCreatorFree::~ParticleCreatorFree()
{
}



//--- Methods ---


ostream &ParticleCreatorFree::write(ostream &s, int shift)
{
	// only ONE of the PCs writes all free particles into ONE .pos file
	if(!toFileDone) {
		s << PutTab(shift) << "<ParticleCreatorFile" << endl;
		s << PutTab(shift+1) << "fileName = \"" << M_SIMULATION->name() << "_restart.pos\"" << endl;
		s << PutTab(shift) << ">" << endl;
		s << PutTab(shift) << "</ParticleCreatorFile>" << endl;

		string fileName(M_SIMULATION->name() + "_restart.pos");
		M_PHASE->writeRestartFile(fileName);
		toFileDone = true;
	}
	
	return s;
}

void ParticleCreatorFree::init()
{
	m_properties.setClassName("ParticleCreatorFree");

  /* Allow unknown properties. Those ones have to be identified later.
	They are used to set the particles degrees of freedom initially. */
	m_properties.allowUnknown();
}


void ParticleCreatorFree::flushParticles()
{
	Phase *phase = M_PHASE;

	if(m_particles.empty())
		throw gError("ParticleCreatorFree::flushParticles", 
		"for " + m_properties.name() + ": the particle list is empty. Is the "
		 "ParticleCreator obsolete? If you REALLY don't think so, make a bug report.");
	size_t counter = 0;
          
	for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++)
	{
		SL_FOR_EACH (Particle, g->second,
			transformVel(*__iSLFE);
			phase->addParticle(*__iSLFE);
			++counter;
		);
	}
	MSG_DEBUG("ParticleCreatorFree::flushParticles for " << m_properties.name()
			, counter << " particles added");
  /* Make sure we don't use too much memory while the simulation
	is running. */
	m_particles.clear();

	MSG_DEBUG("ParticleCreatorFree::flushParticles", "after m_particles.clear()");
/*
	for (size_t i = 0; i < SPACE_DIMS; i++) {
	for (size_t c = 0; c < m_vels[i].size(); c++) {
	delete m_vels[i][c];
	delete m_poss[i][c];
}
}
*/
}

void ParticleCreatorFree::flushParticles(Particle** first_p)
{
	Phase *phase = M_PHASE;
	if(m_particles.empty())
		throw gError
				("ParticleCreatorFree::flushParticles", 
				 "for " + m_properties.name() + ": the particle list is empty. Is the "
						 "ParticleCreator obsolete? If you REALLY don't think so, make a bug report.");
	size_t counter = 0;
          
	for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
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
	     );
	}

	MSG_DEBUG("ParticleCreatorFree::flushParticles for " << m_properties.name()
			, counter << " particles added");

	m_particles.clear();

	MSG_DEBUG("ParticleCreatorFree::flushParticles", "after m_particles.clear()");
}

void ParticleCreatorFree::setup()
{
  ParticleCreator::setup();

  m_temperature_sqrt = sqrt(m_temperature);

}
