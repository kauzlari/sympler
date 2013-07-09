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



#ifndef __THERMOSTAT_WALL_VEL_H
#define __THERMOSTAT_WALL_VEL_H 

#include "thermostat_with_rng.h"
#include "function_fixed.h"


using namespace std;

/* --- ThermostatWallVel --- */

/*!
 * Base class for all wall thermostat. A wall "thermostat" is in charge of
 * resetting the velocities of the frozen particles in each time step.
 * For example, the direction of the velocities can be randomized to 
 * 
 */
class ThermostatWallVel : public ThermostatWithRng
{
protected:
  /*!
   * The species which to thermalize
   */
  string m_species;
  
  /*!
   * The color corresponding to the species
   */
  size_t m_colour;
  
  /*!
   * A bias velocity by which the velocity component is shifted.in the 
   * direction defined by \a m_bias_dir
   */
  FunctionFixed m_bias_velocity;
  
  /*!
   * Determines the orientation of the bias velocity in dependence of 
   */
  point_t m_bias_dir;
    
  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  ThermostatWallVel(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~ThermostatWallVel() {}

  
  /*!
  * setup function
  */
  virtual void setup();
  
  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase* p) = 0;
};

#endif
