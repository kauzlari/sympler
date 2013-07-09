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



#include "thermostat_wall_vel.h"

#include "phase.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


// #define M_SIMULATION ((Simulation*) m_parent)
// #define M_CONTROLLER M_SIMULATION->controller()
// #define M_PHASE M_SIMULATION->phase()
// #define M_MANAGER M_PHASE->manager()


ThermostatWallVel::ThermostatWallVel(Simulation* sim)
  : ThermostatWithRng(sim)
{
}

void ThermostatWallVel::init()
{
  FUNCTIONFIXEDPC(biasVelocity, m_bias_velocity, 
       "The magnitude of the velocity with which the velocity of the"
           " thermalised particle will be biased. You can include time"
           " dependency with the variable 't' and spatial dependency with"
           " the variables 'x', 'y' and 'z'.");

  POINTPC
      (biasDir, m_bias_dir,
       "Vector used for direction of velocity biasing. The magnitude is" 
           " intended to be defined with the attribute 'biasVelocity'."
           " Non-unit vectors given here, will be normalised.");
  
  m_bias_dir.x = m_bias_dir.y = 0;
  m_bias_dir.z = 1;

  m_bias_velocity.addVariable("t");
  m_bias_velocity.addVariable("x");
  m_bias_velocity.addVariable("y");
  m_bias_velocity.addVariable("z");
  m_bias_velocity.setExpression("0");
}

void ThermostatWallVel::setup()
{
  ThermostatWithRng::setup();
  
  m_bias_dir /= m_bias_dir.abs();
}
