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



#include "thermostat_vels.h"

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


const Callable_Register<ThermostatVels> thermostat_vels("ThermostatVels");

/* --- ThermostatVels --- */

ThermostatVels::ThermostatVels(Simulation* sim)
  : Thermostat(sim), m_colour(11111111)
{
  init();
}


void ThermostatVels::init()
{
  m_properties.setClassName("ThermostatVels");

  m_properties.setDescription(
    "Velocity rescaling thermostat."
  );

  DOUBLEPC
      (temperature, temperature, 0,
       "Temperature to thermalize to.");
  temperature = 1;

  STRINGPC
    (species, m_species,
     "Species this thermostat should work on.");

  m_species = "UNDEF";


}


void ThermostatVels::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ThermostatVels::setup", "Attribute 'species' was not defined!");

  m_colour = M_MANAGER->getColour(m_species);
}

void ThermostatVels::thermalize(Phase* phase)
{
    map< size_t, double > velSquare, scale;

    // compute centre of mass velocity
    if(phase -> velCMIsOld()) phase -> computeVelCM();

    // subtract velCM from velocities, because scaling factor can only be computed afterwards

    size_t group;

    FOR_EACH_FREE_PARTICLE_C
      (phase,
       m_colour,
       // velCM per group should be OK, it is questionable whether per colour makes more sense in general
       __iSLFE->v -= phase->centerOfMassVelocity(__iSLFE->g);
       // Compute scaling factor from velSquare
       velSquare[__iSLFE->g] += __iSLFE->v.absSquare();
       group = __iSLFE->g;
      );

    for (map< size_t, double >::iterator g = velSquare.begin(); g != velSquare.end(); g++) {
      g->second /= phase -> returnNofPartC(m_colour);
        scale[g->first] = sqrt(SPACE_DIMS * temperature / g->second);
    }

    // scale to desired temperature
    FOR_EACH_FREE_PARTICLE_C
      (phase,
       m_colour,
       __iSLFE->v *= scale[__iSLFE->g];
       // readd velCM to velocities
       __iSLFE->v += phase->centerOfMassVelocity(__iSLFE->g);
    );

}


