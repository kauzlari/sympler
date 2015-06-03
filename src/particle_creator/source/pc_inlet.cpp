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



#include "pc_inlet.h"

#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"
#include "boundary_with_inlet.h"


/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorInlet> particle_creator_inlet("ParticleCreatorInlet");


#define M_BOUNDARY  ((BoundaryWithInlet *) m_parent)
#define M_PHASE  ((Phase *) M_BOUNDARY->parent())
#define M_MANAGER  M_PHASE->manager()


//---- Constructor/Destructor ----

ParticleCreatorInlet::ParticleCreatorInlet(Boundary *boundary): ParticleCreatorWithRngPCalc(boundary)
{
  init();
}


void ParticleCreatorInlet::init()
{
  m_properties.setClassName("ParticleCreatorInlet");

  m_properties.setDescription(
    "Creates randomly distributed particles within the inlet. Does only work "
    "with BoundarySTL right now."
  );
    
  INTPC(rampSteps, m_nSteps, 0,
        "Number of timesteps for increasing the density linearly from 'initDensity' to 'density'. If rampSteps = 1, the attribute 'initDensity is ignored.'");
  
  DOUBLEPC
      (density, m_density, 0,
       "Desired particle density in the inlet.");
  
  DOUBLEPC
      (initDensity, m_initDensity, 0,
       "This is used as an initial density when a time ramp is used. See also attribute 'rampSteps'");
  
  m_density = 10;
  m_nSteps = 1;
  m_initDensity = HUGE_VAL;
}


void ParticleCreatorInlet::setup()
{
  ParticleCreatorWithRngPCalc::setup();

//   m_temperature_sqrt = sqrt(m_temperature);

  if(m_parent->className() != "BoundaryWithInlet")
    throw gError("ParticleCreatorInlet::setup", "This object can only be used with a Boundary that accepts an inlet");
  
/*  if (typeid(*m_parent) != typeid(BoundarySTL))
    throw gError("ParticleCreatorInlet::setup", "This object can only be used with BoundarySTL.");*/
}


//---- Methods ----


void ParticleCreatorInlet::createParticles()
{    
//   point_t offset, size;
  size_t n_particles, n_initParticles;
  WallContainer &container = M_BOUNDARY->inletBaseSurface();

  m_inlet_length = M_BOUNDARY->inletLength();
  m_inlet_normal = M_BOUNDARY->inletNormal();

  m_inlet_surface = 0;

  list<Wall*> &walls = container.walls();
  for (list<Wall*>::iterator i = walls.begin(); i != walls.end(); ++i) {
/*    MSG_DEBUG
      ("ParticleCreatorInlet::createParticles",
      (*i)->toString());*/

    assert(typeid(*(*i)) == typeid(WallTriangle));

    WallTriangle *wall = (WallTriangle*) *i;
    double surface = wall->side(0).cross(wall->side(1)).abs()/2;
    m_inlet_surface += surface;

    m_walls.push_back(make_pair(wall, surface));
  }
  m_inlet_volume = m_inlet_surface * m_inlet_length;

  MSG_INFO
    ("ParticleCreatorInlet::createParticles",
     "Volume of the inlet is " << m_inlet_volume << ".");
  
  n_particles = (size_t) (m_density*m_inlet_volume);

  initTransform();

  assert(m_createList.empty());

  // do we have a ramp?
  if(m_nSteps > 1)
  {
    if(m_initDensity == HUGE_VAL) 
      throw gError("ParticleCreatorInlet::createParticles", 
                   "Attribute 'initDensity' was not correctly defined!");
    n_initParticles = (size_t) (m_initDensity*m_inlet_volume);
  
    MSG_INFO
        ("ParticleCreatorInlet::createParticles",
         "Creating " << n_particles << " particles during " << m_nSteps << " timestep(s), starting with an initial number of " << n_initParticles << ".");
    
    size_t diff = n_particles - n_initParticles;
    if(diff < m_nSteps) MSG_INFO("ParticleCreatorInlet::createParticles", "Warning: Number of particles needed to increase 'initDensity' to 'density' is lower than 'rampSteps'. Since at least one particle per timestep must be created, the duration of the ramp in timesteps reduces to the number of particles that still have to be created.");
    size_t perStep = diff/m_nSteps; // integer division
    size_t remainder = diff - perStep*m_nSteps;
    // later we read and pop from the front, so the bigger ones must be "pushed_back"
    // the first list entry is skipped because it is done immediately below
    for(size_t i = 1; i < remainder; ++i) m_createList.push_back(perStep+1);
    for(size_t i = remainder; i < m_nSteps; ++i) m_createList.push_back(perStep);
    createMoreParticlesAtStart(n_initParticles+perStep+1);
  }
  else
  {
    MSG_DEBUG
        ("ParticleCreatorInlet::createParticles",
         "Creating " << n_particles << " particles.");
    
    createMoreParticlesAtStart(n_particles);
  }
  /* We create the particle directly */
//  scaleVels();

//  flushParticles();
}


point_t ParticleCreatorInlet::randomPosition()
{
//  RandomNumberGenerator m_rng;
  double r1, r2, h;
  point_t s1, s2, pos;
  double surface = m_rng.uniform()*m_inlet_surface;
  double cs = 0;
  vector< pair<WallTriangle*, double> >::iterator i = m_walls.begin();
  
  cs += i->second;
  while (i != m_walls.end() && cs < surface) {
    ++i;

    assert(i != m_walls.end());

    cs += i->second;
  }

  assert(i != m_walls.end());

  s1 = i->first->side(0);
  s2 = -i->first->side(2);

  r1 = m_rng.uniform();
  r2 = m_rng.uniform();
  h = m_rng.uniform();

  if (r1 + r2 > 1) {
    pos = (1-r1)*s1+(1-r2)*s2;
  } else {
    pos = r1*s1+r2*s2;
  }
//   MSG_DEBUG("ParticleCreatorInlet::randomPosition", "pos before h(=" << h << "): " << pos);
          
  pos += i->first->corner(0) - h*m_inlet_length*m_inlet_normal;
      
//   assert(M_BOUNDARY->isInside(pos));

//   MSG_DEBUG("ParticleCreatorInlet::randomPosition", "m_inlet_length = " << m_inlet_length << ", m_inlet_normal = " << m_inlet_normal);
//   MSG_DEBUG("ParticleCreatorInlet::randomPosition", "pos after h: " << pos);
  
  return pos;
}

// void ParticleCreatorInlet::createMoreParticles()
// {
//   if(!m_createList.empty())
//   {
//     int toCreate = m_createList.front();
//     for(size_t i = 0; i < toCreate; ++i) createParticle(true);
//     m_createList.pop_front();
//   }
// }

void ParticleCreatorInlet::createParticle(bool assign_to_cell)
{
  Cell *c;
  Particle p;
  Particle *new_p;
//   RandomNumberGenerator m_rng;
  
  p.r = randomPosition();
  p.setColour(m_colour);  
  
  transformPos(p);
  
  c = M_MANAGER->findCell(p.r);

  if(!c)
    throw gError("ParticleCreatorInlet::createParticle", "No cell for particle found. Coordinates = ("  + ObjToString(p.r.x) + ", " + ObjToString(p.r.y) + ", " + ObjToString(p.r.z));
  
  p.g = c->group();
  p.dt = 0;
  
  for (int w = 0; w < SPACE_DIMS; w++)
    p.v[w] = m_rng.normal(m_temperature_sqrt);

  transformVel(p);
  
  new_p = M_PHASE->addParticle(p);
/*  MSG_DEBUG("ParticleCreatorInlet::createParticle", "newp: " << new_p->r << endl << new_p->v);*/
  if (assign_to_cell) 
    c->injectFree(m_colour, new_p);
  
/*  Particle &new_p = m_particles[p.g].newEntry();
  new_p.setColour(p.c);
  new_p.r = p.r;
  new_p.v = p.v;
  new_p.g = p.g;
  new_p.tag = p.tag;*/
}


void ParticleCreatorInlet::deleteParticle()
{
  /* Fixme!!! Slow as crap! */
  //RandomNumberGenerator m_rng;
  Cell *c = NULL;
  size_t n, i;
  list<Particle*>::iterator p;

  while (c == NULL) {
    c = m_region->cellAtPos(randomPosition());

    if (c->particles(m_colour).size() == 0) {
      c = NULL;
    } else {
      n = (size_t) 
            (m_rng.uniform() * c->particles(m_colour).size() );

      p = c->particles(m_colour).begin();
      i = 0;
      while (i < n) {
        p++;
        i++;
      }

      c->eraseParticle(p);
    }
  }
}


void ParticleCreatorInlet::scaleVels()
{
  for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
    g->second.scaleVels(m_temperature);
  }
}


