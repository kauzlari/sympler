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


#include "pc_lattice_frozen.h"

#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"


/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorLatticeFrozen> particle_creator_lattice_frozen("ParticleCreatorLatticeFrozen");

#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())

//---- Constructor/Destructor ----

ParticleCreatorLatticeFrozen::ParticleCreatorLatticeFrozen(Boundary *boundary): ParticleCreatorFreePCalc(boundary)
{
    init();
}



//---- Methods ----

/* The ParticleCreatorLatticeFrozen fills the complete space. If one has another PCLF with a smaller box
 * size, this becomes a problem since particles fall out; this gives an Exception. If the size of the
 * old box is smaller than that of the new one, the box size is adjusted */  
void ParticleCreatorLatticeFrozen::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend)
{
  bool lattice_defined = false, dd_defined = false;

  MSG_DEBUG("ParticleCreatorLatticeFrozen::adjustBoxSize", "size = " << size);
    
  if (m_nlattice_points[0] != -1 &&
      m_nlattice_points[1] != -1 &&
      m_nlattice_points[2] != -1) {
    for (int i = 0; i < SPACE_DIMS; i++)
      if (m_nlattice_points[i] < 1)
        throw gError("ParticleCreatorLatticeFrozen::adjustBoxSize", "At least one lattice point "
                     "in each direction is needed.");
        
    MSG_DEBUG("ParticleCreatorLatticeFrozen::adjustBoxSize", "m_nlattice_points = (" 
              << m_nlattice_points[0] << ", " << m_nlattice_points[1] << ", " 
              << m_nlattice_points[2] << ")");
        
    lattice_defined = true;
  }
    
  if (m_distance != HUGE_VAL) {
    if (m_distance <= 0)
	    throw gError("ParticleCreatorLatticeFrozen::adjustBoxSize", "Distance must be bigger than zero.");
        
    MSG_DEBUG("ParticleCreatorLatticeFrozen::adjustBoxSize", "distance = " << m_distance);
        
    dd_defined = true;
  }
    
  if (m_density != HUGE_VAL) {
    if (dd_defined)
	    throw gError("ParticleCreatorLatticeFrozen::adjustBoxSize", "Please specify only one of: " \
                   "distance, density");
        
    MSG_DEBUG("ParticleCreatorLatticeFrozen::adjustBoxSize", "density = " << m_density);
        
    m_distance = pow(1/m_density, (1./3));
    MSG_DEBUG("ParticleCreatorLatticeFrozen::adjustBoxSize", "distance = " << m_distance);        
    dd_defined = true;
  }
    
  if (!lattice_defined && !dd_defined)
	  throw gError("ParticleCreatorLatticeFrozen::adjustBoxSize", "Please specify one of: " \
                 "distance, density, nLatticeP* (with '*' from 'X', 'Y', 'Z')");
					 
    /* If a number of lattice points AND a distance or density AND forceBoundarySize 
    was defined, then the boundary is stretched. */
    if (lattice_defined && dd_defined && m_force_boundary_size) {
      point_t tmpSize;
      bool mustAbort = false;
      for (int i = 0; i < SPACE_DIMS; i++)
      {
        tmpSize[i] = m_nlattice_points[i] * m_distance;
        if(tmpSize[i] < 0.999*size[i])
        {
		MSG_INFO("ParticleCreatorLatticeFrozen::adjustBoxSize", "You are not allowed to"
          " reduce a previously defined box size with settings in this" 
			  " ParticleCreatorLatticeFrozen. This seems to be the case for the " << string(1, 'x'+i) 
          << "-direction:\nold size: " << size[i] << "\nnew size: " << tmpSize[i] 
          << "\nHave to abort soon.");
          mustAbort = true;
        }
      }
      if(mustAbort) throw gError("ParticleCreatorLatticeFrozen::adjustBoxSize: aborting." 
      " Presumably, you can solve this problem either by defining a smaller box size in"
      " the boundary, or by avoiding two ParticleCreatorLattice to force a box size with" 
      " the second forced size being smaller than the first.");
      
      size = tmpSize;
      MSG_INFO("ParticleCreatorLatticeFrozen::adjustBoxSize",
                 "Stretching boundary: size = " << size);
    } 
    else if (!lattice_defined) 
    {
	    MSG_DEBUG("ParticleCreatorLatticeFrozen::adjustBoxSize", "Lattice size not defined. "
          "Will be computed, when box size is known.");
        for (int i = 0; i < SPACE_DIMS; i++) {
/*            int n = (int) (size[i] / m_distance + 0.5);
            if (n < 1)
                n = 1;*/
            m_nlattice_points[i] = /*n*/-1;
        }
    } 
}


void ParticleCreatorLatticeFrozen::createParticles()
{    
	RandomNumberGenerator m_rng;
  	point_t box_size, spacing, offset;
  	double density;
  	size_t seed;
  	/* We need the manager to check which group a particle has. */
  	ManagerCell *manager = M_PHASE->manager();
  	box_size = M_BOUNDARY->boundingBox().size();
  	MSG_DEBUG("ParticleCreatorLatticeFrozen::createParticles", "box_size=" << box_size);
  
  	if (m_nlattice_points[0] == -1 &&
      	m_nlattice_points[1] == -1 &&
      	m_nlattice_points[2] == -1) 
  	{ 
    		for (int i = 0; i < SPACE_DIMS; i++)
		{
        		int n = (int) (box_size[i] / m_distance + 0.5);
        		if (n < 1)
            		n = 1;
        		m_nlattice_points[i] = n;
    		}
  	}
  	offset = M_BOUNDARY->boundingBox().corner1;
  	initTransform();
  	if (m_randomize)
	{
    		seed = 2*getpid() + 1;
  	}
	else
	{
    		seed = 1;
    		MSG_INFO("ParticleCreatorLatticeFrozen::randomizeVels", "randomize = no --> Numbers won't be random.");
  	}

  	density = 1;
  	for (int i = 0; i < SPACE_DIMS; i++) 
	{
    		spacing[i] = m_distance;
    		offset[i] += spacing[i]/2;
    		density *= spacing[i];
  	}    
  	MSG_DEBUG("ParticleCreatorLatticeFrozen::createParticles", "spacing = " << spacing);

  /* The density that is generated may differ from the density wished in the configuration
     file because the particle are equally spaced in space and have the same distance from
     the boundaries (of the cube, of course). That means the REAL density can still differ
     from this value */
  	density = 1/density;
  	MSG_DEBUG("ParticleCreatorLatticeFrozen::createParticles", "density = " << density);
  	MSG_DEBUG("ParticleCreatorLatticeFrozen::createParticles", "lattice = " << m_nlattice_points);

  /* All particles sit in a surrounding box of identical size. This is of course assuming the
     Boundary is exactly as big as returnMaxBoxSize proposes. */

  // creation of FREE particles
  	MSG_DEBUG("ParticleCreatorLatticeFrozen::createParticles", "Creating frozen lattice particles.");
                
  	for (int i = 0; i < m_nlattice_points[0]; i++)
    	for (int j = 0; j < m_nlattice_points[1]; j++)
      	for (int k = 0; k < m_nlattice_points[2]; k++) 
	{
        point_t pos = { { {
                i*spacing[0]+offset[0],
                j*spacing[1]+offset[1],
                k*spacing[2]+offset[2]
            } } };
        Cell *c;
        Particle p;
	
        p.r = pos;
        p.setColour(m_colour);

        transformPos(p);

        c = manager->findCell(p.r);
	
        if (c)
	{
        	if (M_BOUNDARY->isInside(p.r)) 
		{
            		p.g = c->group();
            		for (int w = 0; w < SPACE_DIMS; w++)
              			p.v[w] = m_rng.uniform() - 0.5;
            		m_particles[p.g].newEntry() = p;
          	}
        }
        }
  	scaleVels();
	Phase *phase = M_PHASE;
	if(m_particles.empty())
		throw gError("ParticleCreatorLatticeFrozen::createParticles", 
		": the particle list is empty. Is the "
		 "ParticleCreator obsolete? If you REALLY don't think so, make a bug report.");
	size_t counter = 0;
          
	for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++)
	{
		SL_FOR_EACH (Particle, g->second,
			transformVel(*__iSLFE);
			phase->addFrozenParticle(*__iSLFE);
			++counter;
		);
	}

	MSG_DEBUG("ParticleCreatorLatticeFrozen::flushParticles for " << m_properties.name()
			, counter << " particles added");
  	/* Make sure we don't use too much memory while the simulation is running. */
	m_particles.clear();

	MSG_DEBUG("ParticleCreatorLatticeFrozen::flushParticles", "after m_particles.clear()");
}



void ParticleCreatorLatticeFrozen::init()
{
	m_properties.setClassName("ParticleCreatorLatticeFrozen");

  m_properties.setDescription(
    "Generates frozen particles sitting on a simple cubic lattice. Note that not all"
    " attributes can be set independently. "
    "For example giving a certain density fixes the distance.\n"
    "Note that the size of the simulation box can be modified by this"
    " ParticleCreator. For example, giving "
    "a density and a number of lattice points fixes the size. If this conflicts with the size "
    "given in the Boundary, the Boundary is being resized to reflect the information given here."
  );
    
  DOUBLEPC
    (temperature, m_temperature, 0,
     "Initial temperature of the particles. Note that this "
     "can be overriden by setting vel*.");

  BOOLPC
    (randomize, m_randomize,
     "Set to 'no' to get the same random numbers for every simulation run.");
    
  for (int i = 0; i < SPACE_DIMS; i++) {
    m_properties.addProperty
      ("nLatticeP" + string(1, 'X'+i), PropertyList::INT, &m_nlattice_points[i],
       new PLCIntGreater(0),
       "Number of lattice points (which means particles created) in " + string(1, 'x'+i) + "-"
       "direction.");
    m_nlattice_points[i] = -1;
  }
    
  DOUBLEPC
    (distance, m_distance, 0,
     "Spacing of the lattice points.");
  DOUBLEPC
    (density, m_density, 0,
     "Density of the particles.");

  BOOLPC(forceBoundarySize, m_force_boundary_size,
         "Fill the whole Boundary with particles and the density/lattice information "
         "given here. If necessary, the Boundary is going to be resized. You are only " 
         "allowed to INCREASE a previously defined box size.");
    
  m_randomize = false;
  m_distance = m_density = HUGE_VAL;
  m_temperature = 1;
  m_force_boundary_size = true;
}


void ParticleCreatorLatticeFrozen::setup()
{
  ParticleCreatorFree::setup();

  ParticleCreator::s_createFrozenParts = true;
}


void ParticleCreatorLatticeFrozen::scaleVels()
{
    for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
    //    for (map<int, list<Particle> >::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
        g->second.scaleVels(m_temperature);
    }

}
