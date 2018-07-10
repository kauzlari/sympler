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



#include "f_wall_repulsion.h"

#include "random.h"
#include "simulation.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()


const GenFTypeConcr<FWallRepulsion> f_wall_repulsion("FWallRepulsion");

//---- Constructors/Destructor ----

FWallRepulsion::FWallRepulsion(Simulation *simulation)
  : FParticle(simulation)
{
    init();
}


FWallRepulsion::~FWallRepulsion()
{
}


void FWallRepulsion::init()
{
  m_properties.setClassName("FWallRepulsion");

  m_properties.setDescription
    ("Heat flux from two parallel walls to the fluid."
     );

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Cut-off distance from the wall.");

  DOUBLEPC
    (force, m_force, 0,
     "Heat conduction coefficient.");

  INTPC
    (wallDir, m_wall_dir, -1,
     "Direction in which to find the wall: 0 = x, 1 = y, 2 = z.");

  m_properties.addProperty
    ("leftWall", PropertyList::DOUBLE, &m_left_wall, NULL,
     "Position of the wall to the left.");

  m_properties.addProperty
    ("rightWall", PropertyList::DOUBLE, &m_right_wall, NULL,
     "Position of the wall to the right.");

  m_wall_dir = 0;

  m_cutoff = 1;
  m_force = 1;

  m_left_wall = -10;
  m_right_wall = 10;

  m_is_pair_force = false;
  m_is_particle_force = true;
}


//---- Methods ----

void FWallRepulsion::computeForces(Particle* part, int force_index)
{
  Phase *phase = M_PHASE;
  point_t n = {{{ 0, 0, 0 }}};

  n[m_wall_dir] = 1;

  FOR_EACH_PARTICLE_C
    (phase,
     m_colour,

     double weight = 0;

     if (__iSLFE->r[m_wall_dir] < m_left_wall+m_cutoff) {
       weight = 1-(__iSLFE->r[m_wall_dir]-m_left_wall)*m_rcinv;
     } else if (__iSLFE->r[m_wall_dir] > m_right_wall-m_cutoff) {
       weight = -(1-(m_right_wall-__iSLFE->r[m_wall_dir])*m_rcinv);
     }

     __iSLFE->force[force_index] += weight * m_force * n;
     );
}


#ifndef _OPENMP
void FWallRepulsion::computeForces(Pairdist* pair, int force_index)
#else
void FWallRepulsion::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FWallRepulsion::computeForces", "Fatal error: do not call FWallRepulsion::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FWallRepulsion::computeForces(int force_index)
{
  throw gError("FWallRepulsion::computeForces", "Fatal error: do not call FWallRepulsion::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FWallRepulsion::setup()
{
  FParticle::setup();

  m_rcinv = 1/m_cutoff;
}


