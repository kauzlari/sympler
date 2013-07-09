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



#include "thermostat_energy_rescaling.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

const Callable_Register<ThermostatEnergyRescaling> thermostat_energy_rescaling("ThermostatEnergyRescaling");


ThermostatEnergyRescaling::ThermostatEnergyRescaling(Simulation* sim)
  : Thermostat(sim)
{
  init();
}


void ThermostatEnergyRescaling::init()
{
  m_properties.setClassName("ThermostatEnergyRescaling");

  m_properties.setDescription
    ("Rescale to get a predefined average energy.");

  STRINGPC
    (species, m_species,
     "Species this thermostat should work on.");

  DOUBLEPC
    (internalEnergy, m_internal_energy, 0,
     "Reset average internal energy to this value.");

  m_species = "UNDEF";
  m_internal_energy = -1;
}


void ThermostatEnergyRescaling::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ThermostatEnergyRescaling::setup", "Attribute 'species' was not defined!");

  m_colour = M_MANAGER->getColour(m_species);

  m_ie = 
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species);

  if (!m_ie)
    throw gError
      ("TheromstatEnergyRescaling::read",
       "You cannot use this object without IntegratorEnergy for the"
           " corresponding species.");

  m_energy_offset =
    Particle::s_tag_format[m_colour].attrByName("internal_energy").offset;
}


void ThermostatEnergyRescaling::thermalize(Phase* phase)
{
  double e = 0;
  size_t n = 0;
  double factor;

  FOR_EACH_PARTICLE_C
    (phase,
     m_colour,
     e += __iSLFE->tag.doubleByOffset(m_energy_offset);
     n++;
     );

  e /= n;

  factor = m_internal_energy / e;

  MSG_DEBUG
    ("ThermostatEnergyRescaling::thermalize",
     "factor = " << factor << ", <e> = " << m_internal_energy);

  FOR_EACH_PARTICLE_C
    (phase,
     m_colour,
     __iSLFE->tag.doubleByOffset(m_energy_offset) *= factor;
     );
}
