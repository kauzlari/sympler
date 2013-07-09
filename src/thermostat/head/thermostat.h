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


#ifndef __THERMOSTAT_H
#define __THERMOSTAT_H 

#include "random.h"
#include "general.h"
#include "callable.h"
#include "colour_pair.h"

using namespace std;

class Phase;
class Simulation;

/*!
 * A \a Thermostat is the base class for each thermostat, i.e., modules that keep the
 * temperature of the system constant.
 */
class Thermostat : public Callable
{
 protected:
  /*!
   * Time step when to activate this thermostat
   */
  int m_activate_at;

  /*!
   * Time step when to deactivate this thermostat
   */
  int m_deactivate_at;

  /*!
   * The interval when to call the thermostats
   */
  size_t m_interval;

  /*!
   * The current time step
   */
  size_t m_step;

  /*!
   * Is the thermostat active?
   */
  bool m_active;

  /*!
   * Initialize the property list
   */
  void init();

 public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  Thermostat(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~Thermostat();

  /*!
   * This function decides when to call the \a thermalize method,
   * i.e., only at certain intervals given by \a m_interval, and
   * only if the thermostat is active.
   */
  virtual void call(size_t timestep);

  /*!
   * The routing thermalizing the system
   */
  virtual void thermalize(Phase* p) = 0;

  /*!
   * Set \a m_active to false
   */
  virtual void setup();
};


#endif

