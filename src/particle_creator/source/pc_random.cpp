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



#include "pc_random.h"

#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"


/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorRandom> particle_creator_random("ParticleCreatorRandom");


#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())



//---- Constructor/Destructor ----

ParticleCreatorRandom::ParticleCreatorRandom(Boundary *boundary): ParticleCreatorWithRngPCalc(boundary)
{
  init();
}



//---- Methods ----


void ParticleCreatorRandom::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend)
{
  bool n_defined = m_nparticles > 0, d_defined = m_density > 0;

  MSG_DEBUG("ParticleCreatorRandom::adjustBoxSize", "Old size = " << size);

  if (!n_defined && !d_defined)
    throw gError
      ("ParticleCreatorRandom::adjustBoxSize",
       "Please define either nParticles or density.");

  if (n_defined && d_defined) {
    double stretchFactor = pow(((m_nparticles/m_density) / M_BOUNDARY->cuboidVolume()), (1./3));
    for(size_t i = 0; i < SPACE_DIMS; ++i)
      size[i] *= stretchFactor;
//     size.assign(pow(m_nparticles/m_density, (1./3)));
  } else {
    if (n_defined)
      m_density = m_nparticles/M_BOUNDARY->cuboidVolume();
    else
      m_nparticles = (int) (m_density*M_BOUNDARY->cuboidVolume() + 0.5);
  }

  m_offset = M_BOUNDARY->boundingBox().corner1;
  m_size = size;

  MSG_DEBUG("ParticleCreatorRandom::adjustBoxSize", "New size = " << size);
}


void ParticleCreatorRandom::createParticles()
{    
//   int seed;
  ManagerCell *manager;
  double m_sqrt_t = sqrt(m_temperature);
  
  manager = M_PHASE->manager();

  MSG_DEBUG
    ("ParticleCreatorRandom::createParticles",
     "m_nparticles = " << m_nparticles << ", m_density = " << m_density);

  initTransform();

  int i = 0;
  
  while(i < m_nparticles)
  
 /* for (int i = 0; i < m_nparticles; i++)*/ {
    
    // ++i is OK here, because the PC creates particles for the whole bounding box. So 
    // it is OK if a particle, which is outside the geometry will not be created
    ++i;
    
    point_t pos = { { {
            m_rng.uniform()*m_size.x+m_offset.x,
            m_rng.uniform()*m_size.y+m_offset.y,
            m_rng.uniform()*m_size.z+m_offset.z
        } } };

    Particle p;
    // next is necessary before transformPos!
    p.setColour(m_colour);  
    p.r = pos;
    transformPos(p);
      
    if(M_BOUNDARY->isInside(p.r))
    {
        /*MSG_DEBUG("ParticleCreatorRandom::createParticles", "INSIDE no " << i << ", goal = " << m_nparticles << ", p = " << pos);*/
      
      Cell *c;
      c = manager->findCell(p.r);
  
      if (c) {
        p.g = c->group();
  
	//        for (int w = 0; w < SPACE_DIMS; w++)
	//          p.v[w] = m_rng.uniform() - 0.5;

        for (int w = 0; w < SPACE_DIMS; w++)
          p.v[w] = m_sqrt_t * (2*m_rng.uniform()-1);
  
        Particle &new_p = m_particles[p.g].newEntry();
        new_p.setColour(p.c);
        new_p.r = p.r;
        new_p.v = p.v;
        new_p.g = p.g;
        new_p.tag = p.tag;
//         MSG_DEBUG("ParticleCreatorRandom::createParticles", "new_p.v = " << new_p.v);
      }
      else
      {
        cout << "The following error occured:" << endl;
        cout << "[ParticleCreatorRandom::createParticles] no cell for particle with position " << p.r << endl;
        abort();
      }
    }
//     else MSG_DEBUG("ParticleCreatorRandom::createParticles", "NOT INSIDE, p = " << pos);
  }

  scaleVels();

  flushParticles();
}



void ParticleCreatorRandom::init()
{
  m_properties.setClassName("ParticleCreatorRandom");

  m_properties.setDescription(
    "Creates randomly distributed particles. Does only work with a cuboid volume!"
  );
    
  DOUBLEPC
    (temperature, m_temperature, 0,
     "Initial temperature of the particles. Note that this "
     "can be overriden by setting vel*.");

  INTPC
    (nParticles, m_nparticles, -2,
     "Number of particles.");
  
  DOUBLEPC
    (density, m_density, -2,
     "Density of the particles.");
  
  m_density = -1;
  m_nparticles = -1;
  m_temperature = 1;
}


void ParticleCreatorRandom::scaleVels()
{
  for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
    g->second.scaleVels(m_temperature);
  }
}
