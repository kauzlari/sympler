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


#include "thermostat_wall_vel_iso.h"

#include "phase.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<ThermostatWallVelIso> 
  thermostat_wall_vel_iso("ThermostatWallVelIso");


ThermostatWallVelIso::ThermostatWallVelIso(Simulation* sim)
  : ThermostatWallVel(sim)
{
  init();
}


void ThermostatWallVelIso::init()
{
  ThermostatWallVel::init();
  
  m_properties.setClassName("ThermostatWallVelIso");

  m_properties.setDescription(
    "Sets the velocity of wall particles according to the defined wall temperature.");

  STRINGPC
    (species, m_species,
     "Species this thermostat should work on.");
  m_species = "UNDEF";
     
  DOUBLEPC
    (kBToverM, m_temperature, 0,
     "k_BT/m = <v^2> to thermalize wall-particles to.");
  m_temperature = 1;
}


void ThermostatWallVelIso::setup()
{
  ThermostatWallVel::setup();

  m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);

}


void ThermostatWallVelIso::thermalize(Phase* phase)
{
//  RandomNumberGenerator m_rng;
  double time = M_CONTROLLER->time();
  
  FOR_EACH_FROZEN_PARTICLE
    (phase, m_colour,
     double v = sqrt(SPACE_DIMS*m_temperature);
     double phi = m_rng.uniform()*M_PI;
     double theta = m_rng.uniform()*2*M_PI;
     double sin_phi = sin(phi);

     __iSLFE->v.x = v*sin_phi*cos(theta);
     __iSLFE->v.y = v*sin_phi*sin(theta);
     __iSLFE->v.z = v*cos(phi);
    
     // add the bias
     __iSLFE->v += m_bias_velocity(time, __iSLFE->r.x, __iSLFE->r.y, __iSLFE->r.z)*m_bias_dir;
    );
}
