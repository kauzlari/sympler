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



#ifndef __ENERGY_PUMP_H
#define __ENERGY_PUMP_H 

#include "callable.h"
#include "integrator_energy.h"

using namespace std;

/* --- EnergyPump --- */


/*!
 * Pump a certain amount of energy into the internal energy of
 * each particle.
 */
class EnergyPump : public Callable
{
protected:
  /*!
   * Species to pump energy into
   */
  string m_species;

  /*!
   * Color to pump energy into
   */
  size_t m_colour;

  /*!
   * Tag offset of the internal energy
   */
  size_t m_energy_offset;

  /*!
   * Current simulation step
   */
  size_t m_step;

  /*!
   * Energy to pump into each particle in the interval \a m_interval
   */
  double m_energy_per_particle;

  /*!
   * Pump energy every \a m_interval time steps
   */
  size_t m_interval;

  /*!
   * Initialize property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param sim Pointer to the simulation object
   */
  EnergyPump(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~EnergyPump() {}

  /*!
   * Increases the energy for each particle
   */
  virtual void call(size_t timestep);

  /*!
   * Get the color and energy offset
   */
  virtual void setup();
};

#endif
