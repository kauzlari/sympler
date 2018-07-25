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



#include "f_coriolis.h"
#include "manager_cell.h"
#include "simulation.h"

#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<FCoriolis> f_coriolis("FCoriolis");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

FCoriolis::FCoriolis(Simulation *simulation): FParticle(simulation)
{
    init();
}


FCoriolis::~FCoriolis()
{
}



//---- Methods ----


void FCoriolis::computeForces(Particle* part, int force_index)
{
  Phase *phase = ((Simulation *) m_parent)->phase();
  double freq = m_frequency(M_CONTROLLER->time());

  MSG_DEBUG("FCentrigual::computeForces", "freq = " << freq);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     __iSLFE->force[force_index] += freq*m_rotation_axis.cross(__iSLFE->v);
    );
}


#ifndef _OPENMP
void FCoriolis::computeForces(Pairdist* pair, int force_index)
#else
void FCoriolis::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
    throw gError("FCoriolis::computeForces", "Fatal error: do not call FCoriolis::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FCoriolis::computeForces(int force_index)
{
    throw gError("FCoriolis::computeForces", "Fatal error: do not call FCoriolis::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FCoriolis::init()
{
  m_properties.setName("FCoriolis");

  m_properties.setDescription("Centrifugal force.");

  POINTPC
    (rotationAxis, m_rotation_axis,
     "Axis of rotation. Note: This value is being normalized! Use frequency instead.");

  FUNCTIONFIXEDPC
    (frequency, m_frequency,
     "Frequency of the rotation.");

  m_rotation_axis.x = m_rotation_axis.y = 0;
  m_rotation_axis.z = 1;

  m_frequency.addVariable("t");

  m_is_pair_force = false;
  m_is_particle_force = false;
}


void FCoriolis::setup()
{
  FParticle::setup();

  m_rotation_axis /= m_rotation_axis.abs();
}

