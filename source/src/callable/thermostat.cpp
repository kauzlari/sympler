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



#include "thermostat.h"

#include "phase.h"
#include "random.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


/* --- Thermostat --- */

Thermostat::Thermostat(Simulation* sim)
  : Callable(sim)
{
  init();
}


Thermostat::~Thermostat()
{
}


void Thermostat::init()
{
  m_properties.setClassName("Thermostat");

  INTPC
    (activateAt, m_activate_at, -1,
     "Timestep at which the thermostat should be activated.");

  INTPC
    (deactivateAt, m_deactivate_at, -2,
     "Timestep at which the thermostat should be deactivated. -1 means never deactivate.");

  INTPC
    (interval, m_interval, 0,
     "Call interval for the termostat.");

  m_activate_at = 0;
  m_deactivate_at = -1;
  m_interval = 1;
}


void Thermostat::setup()
{
  NodeManyChildren::setup();

  m_active = false;
}


void Thermostat::call(size_t timestep)
{
  /*  MSG_DEBUG
    ("Thermostat::call",
    m_properties.name() + " active = " << m_active);*/

  if (m_active) {
    if ((int) timestep == m_deactivate_at) {
      m_active = false;
    } else {
      m_step--;

      if (!m_step) {
	thermalize(M_PHASE);
	m_step = m_interval;
      }
    }
  } else if ((int) timestep == m_activate_at) {
    m_active = true;
    m_step = m_interval;
    thermalize(M_PHASE);
  }
}

