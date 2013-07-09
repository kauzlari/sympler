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



#include "energy_pump.h"


#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<EnergyPump> energy_pump("EnergyPump");


EnergyPump::EnergyPump(Simulation* sim)
  : Callable(sim)
{
  init();
}


void EnergyPump::init()
{
  m_properties.setClassName("EnergyPump");

  m_properties.setDescription(
    "Pumps a certain amount of energy into each particles internal energy."
  );

  STRINGPC
    (species, m_species,
     "Species this pump should work on.");

  DOUBLEPC
    (energyPerParticle, m_energy_per_particle, 0,
     "Energy to pump per particle.");

  INTPC
    (interval, m_interval, 0,
     "Pump every <interval> timesteps.");

  m_species = "fluid";
  m_energy_per_particle = 1;
  m_interval = 100;
}


void EnergyPump::setup()
{
  Callable::setup();

  m_colour = M_MANAGER->getColour(m_species);

  m_energy_offset = Particle::s_tag_format[m_colour].attrByName("internal_energy").offset;

  m_step = 0;
}


void EnergyPump::call(size_t timestep)
{
  Phase *phase = M_PHASE;

  m_step++;

  if (m_step == m_interval) {
    MSG_DEBUG("EnergyPump::thermalize", "Pumping energy.");

    FOR_EACH_PARTICLE_C
      (phase,
       m_colour,
       __iSLFE->tag.doubleByOffset(m_energy_offset) += m_energy_per_particle;
      );

    m_step = 0;
  }
}


