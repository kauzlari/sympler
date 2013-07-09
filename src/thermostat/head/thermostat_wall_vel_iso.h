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


#ifndef __THERMOSTAT_WALL_VEL_ISO_H
#define __THERMOSTAT_WALL_VEL__ISO_H 

#include "thermostat_wall_vel.h"

using namespace std;

/* --- ThermostatWallVelIso --- */

/*!
 * Randomize the velocities of frozen particles to a pre-defined temperature
 */
class ThermostatWallVelIso : public ThermostatWallVel
{
protected:
  /*!
   * The temperature to which to thermalize to
   */
  double m_temperature;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  ThermostatWallVelIso(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~ThermostatWallVelIso() {}

  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase* p);

  /*!
   * Get the color
   */
  virtual void setup();
};

#endif
