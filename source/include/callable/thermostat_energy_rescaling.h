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



#ifndef __THERMOSTAT_ENERGY_RESCALING_H
#define __THERMOSTAT_ENERGY_RESCALING_H 

#include "thermostat.h"
#include "integrator_energy.h"

using namespace std;

/* --- ThermostatEnergyRescaling --- */


/*!
 * Rescales the internal energies of the particle
 * to match a certain average energy.
 */
class ThermostatEnergyRescaling : public Thermostat
{
protected:
  /*!
   * The species for which to rescale the energy
   */
  string m_species;

  /*!
   * The color of the species
   */
  size_t m_colour;

  /*!
   * Tag offset of the internal energy
   */
  size_t m_energy_offset;
  
  /*!
   * Integrator for the internal energy degree of freedom
   */
  IntegratorEnergy *m_ie;

  /*!
   * Energy to rescale to
   */
  double m_internal_energy;

  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Find the integrator and the tag offset of the internal energy
   */
  virtual void setup();

public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  ThermostatEnergyRescaling(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~ThermostatEnergyRescaling() {}

  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase* p);
};

#endif
