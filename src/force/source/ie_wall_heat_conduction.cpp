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



#include "ie_wall_heat_conduction.h"

#include "random.h"
#include "simulation.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()


const GenFTypeConcr<IEWallHeatConduction> ie_wall_heat_conduction("IEWallHeatConduction");

//---- Constructors/Destructor ----

IEWallHeatConduction::IEWallHeatConduction(Simulation *simulation)
  : FParticle(simulation)
{
    init();
}


IEWallHeatConduction::~IEWallHeatConduction()
{
}


void IEWallHeatConduction::init()
{
  m_properties.setClassName("IEWallHeatConduction");

  m_properties.setDescription
    ("Heat flux from two parallel walls to the fluid."
     );

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Cut-off distance from the wall.");

  DOUBLEPC
    (kappa, m_kappa, 0,
     "Heat conduction coefficient.");

  INTPC
    (wallDir, m_wall_dir, -1,
     "Direction in which to find the wall: 0 = x, 1 = y, 2 = z.");

  m_properties.addProperty
    ("leftWall", PropertyList::DOUBLE, &m_left_wall, NULL,
     "Position of the wall to the left.");

  DOUBLEPC
    (leftTemperature, m_left_temperature, 0,
     "Temperature of the wall to the left.");

  m_properties.addProperty
    ("rightWall", PropertyList::DOUBLE, &m_right_wall, NULL,
     "Position of the wall to the right.");

  DOUBLEPC
    (rightTemperature, m_right_temperature, 0,
     "Temperature of the wall to the right.");

  m_properties.addProperty
    ("spotSize", PropertyList::DOUBLE, &m_spot_size, NULL,
     "Size of the spot to heat.");

  m_wall_dir = 0;

  m_spot_size = -1;

  m_cutoff = 1;
  m_kappa = 1;

  m_left_wall = -10;
  m_left_temperature = 1;
  m_right_wall = 10;
  m_right_temperature = 1;

  m_is_pair_force = false;
  m_is_particle_force = true;
}


//---- Methods ----

void IEWallHeatConduction::computeForces(Particle* part, int force_index)
{
//  RandomNumberGenerator m_rng;
  Phase *phase = M_PHASE;

  if (m_spot_size < 0) {
    FOR_EACH_PARTICLE_C
      (phase,
       m_colour,

       double weight = -1;
       double Tinv;

       if (__iSLFE->r[m_wall_dir] < m_left_wall+m_cutoff) {
	 weight = 1-(__iSLFE->r[m_wall_dir]-m_left_wall)*m_rcinv;
	 Tinv = m_lTinv;
       } else if (__iSLFE->r[m_wall_dir] > m_right_wall-m_cutoff) {
	 weight = 1-(m_right_wall-__iSLFE->r[m_wall_dir])*m_rcinv;
	 Tinv = m_rTinv;
       }

       if (weight > 0) {
	 double f = weight * weight * m_kappa
	   * (m_ie->reciprocalTemperature(*__iSLFE) - Tinv);
	 /* fixme!!! commented out heat fluctuation.
	    because of the large mass of the wall? */
	 //         + weight * m_alpha * m_rng.normal(1) * m_r_sqrt_dt;

	 __iSLFE->tag.doubleByOffset(m_eforce_offset[force_index]) += f;
       }
       );
  } else {
    FOR_EACH_PARTICLE_C
      (phase,
       m_colour,

       double weight = -1;
       double Tinv;

       if (__iSLFE->r[m_wall_dir] < m_left_wall+m_cutoff) {
	 weight = 1-(__iSLFE->r[m_wall_dir]-m_left_wall)*m_rcinv;
	 Tinv = m_lTinv;
       } else if (__iSLFE->r[m_wall_dir] > m_right_wall-m_cutoff) {
	 weight = 1-(m_right_wall-__iSLFE->r[m_wall_dir])*m_rcinv;
	 Tinv = m_rTinv;
       }

       if (weight > 0 &&
	   fabs(__iSLFE->r[m_perp_dir1]) < m_spot_size &&
	   fabs(__iSLFE->r[m_perp_dir2]) < m_spot_size) {
	 double f = weight * weight * m_kappa
	   * (m_ie->reciprocalTemperature(*__iSLFE) - Tinv);
	 /* fixme!!! commented out heat fluctuation.
	    because of the large mass of the wall? */
	 //         + weight * m_alpha * m_rng.normal(1) * m_r_sqrt_dt;

	 __iSLFE->tag.doubleByOffset(m_eforce_offset[force_index]) += f;
       }
       );
  }
}


#ifndef _OPENMP
void IEWallHeatConduction::computeForces(Pairdist* pair, int force_index)
#else
void IEWallHeatConduction::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("IEWallHeatConduction::computeForces", "Fatal error: do not call IEWallHeatConduction::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void IEWallHeatConduction::computeForces(int force_index)
{
  throw gError("IEWallHeatConduction::computeForces", "Fatal error: do not call IEWallHeatConduction::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void IEWallHeatConduction::setup()
{
  FParticle::setup();

  m_rcinv = 1/m_cutoff;

  m_lTinv = 1/m_left_temperature;
  m_rTinv = 1/m_right_temperature;

  m_alpha = sqrt(2*m_kappa);

  m_r_sqrt_dt = 1/sqrt(M_CONTROLLER->dt());

  m_colour = M_MANAGER->getColour(m_species);

  m_perp_dir1 = (m_wall_dir+1)%SPACE_DIMS;
  m_perp_dir2 = (m_wall_dir+2)%SPACE_DIMS;

  /* Get integrators */
  m_ie =
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species);

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
    m_eforce_offset[i] =
        Particle::s_tag_format[m_colour].attrByName("force_internal_energy_" + ObjToString(i)).offset;

  /* fixme!!! Do a real check */
  assert(m_ie);
}


