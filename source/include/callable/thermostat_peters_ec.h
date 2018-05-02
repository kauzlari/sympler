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



#ifndef __THERMOSTAT_PETERS_EC_H
#define __THERMOSTAT_PETERS_EC_H 

#include "thermostat_peters.h"
#include "integrator_energy.h"
#include "weighting_function.h"


using namespace std;

/*!
 * Implementation of the modified energy conserving Peters thermostat.
  * See: Phys. Rev. E 73 037701
 * For the original see \a ThermostatPetersIso and: E. A. J. F. Peters, 
  * Europhys. Lett. 66, 311-317 (2004)
 */

class ThermostatPetersEnergyConserving : public ThermostatPeters
{
protected:
  
  /*!
  * memory offset in the particles, where the energy-DOF will be saved
  */
  pair<size_t, size_t> m_energy_offset;
  
  /*!
  * The two Integrators for the energy-DOF of the two species in \a m_species 
  */
  pair<IntegratorEnergy*, IntegratorEnergy*> m_ie;

  /*!
   * Initialize the property list
  */
  void init();

public:
  
  /*!
  * Constructor
  */
  ThermostatPetersEnergyConserving(Simulation* sim);
  
  /*!
   * Destructor
   */
  virtual ~ThermostatPetersEnergyConserving() {}

  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase* p);
  
  /*!
   * Initialize variables
   */
  virtual void setup();
};

#endif
