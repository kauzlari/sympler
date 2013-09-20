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



#include "boundary.h"
#include "simulation.h"
#include "manager_cell.h"
#include "phase.h"
#include "reflector.h"
#include "particle_creator.h"

// #include "valgrind/memcheck.h"

using namespace std;


#define M_PHASE  ((Phase*) m_parent)
#define M_SIMULATION ((Simulation*) M_PHASE->parent())

REGISTER_SMART_ENUM
(Boundary_Factory,
 "As the name suggests, Boundaries define the boundaries of a simulation region. "
 "They are children of a phase."
);


//---- Constructors/Destructor ----

Boundary::Boundary()
    : NodeManyChildren()/*, m_particle_creator(NULL)*/ 
{
    throw gError("Boundary::Boundary: Shouldn't be called! FATAL!");
    // This should never be called!
}


Boundary::Boundary(Phase *phase)
  : NodeManyChildren((Node*) phase)
{
  for(size_t i = 0; i < SPACE_DIMS ; ++i)
    {
      m_frontFrame[i] = false;
      m_endFrame[i] = false;
      m_proposedSize[i] = -1;		
    }
  m_thickness = -1;
}


Boundary::~Boundary()
{

}



//---- Methods ----

void Boundary::createParticles()
{
  for(vector<ParticleCreator*>::iterator pcIter = m_pcList.begin();
      pcIter != m_pcList.end(); ++pcIter)
    (*pcIter)->createParticles();
}

void Boundary::createMoreParticles()
{
  //  MSG_DEBUG("Boundary::createMoreParticles", "START");
  for(vector<ParticleCreator*>::iterator pcIter = m_pcList.begin();
      pcIter != m_pcList.end(); ++pcIter)
    {
      //      MSG_DEBUG("Boundary::createMoreParticles", "CALLING PC");
      (*pcIter)->createMoreParticles();
    }
}


bool Boundary::isInside(point_t point)
{
    return boundingBox().isInside(point);
}


bool Boundary::isInside(cuboid_t cuboid, const double& range)
{
  return isInside(cuboid.corner1) || isInside(cuboid.corner2);
}


Node *Boundary::instantiateChild(const string &name)
{
    Node* node;

    if (ParticleCreator_Factory::exists(name)) {
	    ParticleCreator* pc;

      node = pc = ParticleCreator_Factory::byName(name).instantiate(this);
	
      m_pcList.push_back(pc);
    } else if (Reflector_Factory::exists(name)) {
      Reflector *reflector;

      node = reflector = Reflector_Factory::byName(name).instantiate(/*NULL*/this);

      m_reflectors.push_back(reflector);
    } else
      throw gError
        ("Boundary::instantiateChild", 
         "Object " + name + " not found in database.");

    return node;
}


void Boundary::setup()
{
  NodeManyChildren::setup();

  if(m_pcList.empty())
    throw gError("Boundary::setup", "No ParticleCreator defined.");

  if(m_reflectors.empty())
    throw gError("Boundary::setup", "No Reflector defined.");

  if (!M_PHASE->pairCreator()) {
    throw gError("Boundary::setup", "No PairCreator defined.");

  m_thickness = M_PHASE->pairCreator()->interactionCutoff();
//   m_thickness = M_SIMULATION->maxCutoff;
}


//---- Cell subdivision ----

void Boundary::setup(Simulation* sim, ManagerCell *mgr)
{
  for(vector<ParticleCreator*>::iterator pcIter = m_pcList.begin();
      pcIter != m_pcList.end(); ++pcIter) { 
    (*pcIter) -> adjustBoxSize(m_proposedSize, m_frontFrame, m_endFrame);
  }
}



Reflector *Boundary::findReflector(const string &name)
{
  Reflector *r = NULL;

  FOR_EACH
    (vector<Reflector*>,
     m_reflectors,
     if ((*__iFE)->name() == name)
       r = *__iFE;
    );
// if(r) MSG_DEBUG("Boundary::findReflector", "found Reflector " + r->name());
  if (!r)
    throw gError
      ("Boundary::findReflector", "Unknown reflector '" + name + "'.");

  return r;
}
