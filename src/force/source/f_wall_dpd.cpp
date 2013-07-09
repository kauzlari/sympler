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



#include "f_wall_dpd.h"

#include "random.h"
#include "simulation.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()


const GenFTypeConcr<FWallDPD> f_wall_dpd("FWallDPD");

//---- Constructors/Destructor ----

FWallDPD::FWallDPD(Simulation *simulation): FParticle(simulation)
{
    init();
}


FWallDPD::~FWallDPD()
{
}


void FWallDPD::init()
{
  m_properties.setClassName("FWallDPD");

  m_properties.setDescription
    ("Models the wall as a large flat DPD particle."
     );

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Cut-off distance from the wall.");

  DOUBLEPC
    (dissipation, m_dissipation, 0,
     "Dissipation constant for the wall-fluid interaction.");

  DOUBLEPC
    (density, m_density, 0,
     "Density of the wall in terms of fluid density.");

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

  m_dissipation = 1;
  m_density = 1;

  m_left_wall = -10;
  m_right_wall = 10;

  m_is_pair_force = false;
  m_is_particle_force = true;
}


//---- Methods ----

void FWallDPD::computeForces(Particle* part, int force_index)
{
  Phase *phase = M_PHASE;
  int other1, other2;

  other1 = 0;
  while (other1 == m_wall_dir)
    other1++;

  other2 = 0;
  while (other2 == m_wall_dir || other2 == other1)
    other2++;

  //  cout << m_wall_dir << " " << other1 << " " << other2 << endl;

  FOR_EACH_PARTICLE_C
    (phase,
     m_colour,

     double zp;
     bool near_wall = false;

     if (__iSLFE->r[m_wall_dir] < m_left_wall+m_cutoff) {
       zp = (__iSLFE->r[m_wall_dir]-m_left_wall)/m_cutoff;
       near_wall = true;
     } else if (__iSLFE->r[m_wall_dir] > m_right_wall-m_cutoff) {
       zp = (m_right_wall-__iSLFE->r[m_wall_dir])/m_cutoff;
       near_wall = true;
     }

     if (near_wall) {
       point_t f;
       double perp;
       double par;

       perp = pow(1-zp,4)*(1+4*zp+10*zp*zp);
       par = pow(1-zp,5)*(1+2*zp);

       f[m_wall_dir] = m_factor*perp*__iSLFE->v[m_wall_dir];
       f[other1] = m_factor*par*__iSLFE->v[other1];
       f[other2] = m_factor*par*__iSLFE->v[other2];

       __iSLFE->force[force_index] += f;

     }
     );
}


#ifndef _OPENMP
void FWallDPD::computeForces(Pairdist* pair, int force_index)
#else
void FWallDPD::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FWallDPD::computeForces", "Fatal error: do not call FWallDPD::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FWallDPD::computeForces(int force_index)
{
  throw gError("FWallDPD::computeForces", "Fatal error: do not call FWallDPD::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FWallDPD::setup()
{
  FParticle::setup();

  m_factor = -m_dissipation*m_density*m_cutoff*m_cutoff/12;
}


