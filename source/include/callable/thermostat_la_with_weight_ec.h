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



#ifndef __THERMOSTAT_LA_WITH_WEIGHT_EC_H
#define __THERMOSTAT_LA_WITH_WEIGHT_EC_H 

#include "thermostat_with_rng.h"
#include "integrator_energy.h"
#include "weighting_function.h"

using namespace std;

/* --- ThermostatLAWithWeightEC --- */

/*!
 * Implementation of the Lowe-Andersen thermostat 
 * (C. P. Lowe, Europhys. Lett. 47, 145-151 (1999)) with energy conservation.
 * This version of the thermostat uses a user defined weighting function.
 */
class ThermostatLAWithWeightEC : public ThermostatWithRng
{
protected:
  /*!
   * The species to thermalize
   */
  pair<string, string> m_species;

  /*!
   * The pair of colors to thermalize
   */
  ColourPair *m_cp;

  /*!
   * The dissipation constant
   */
  double m_dissipation;

  /*!
   * Thermalization probability = dissipation * time step
   */
  double m_probability;

  /*!
   * The cut-off radius, taken from the weighting function
   */
  double m_cutoff;

  /*!
   * Tag offset of the internal energy for the two species
   */
  pair<int, int> m_energy_offset;

  /*!
   * Energy integrators for the two species
   */
  pair<IntegratorEnergy*, IntegratorEnergy*> m_ie;

  /*!
   * The weighting function to be used
   */
  string m_weighting_function;

  /*!
   * Pointer to the implementation of the weighting function
   */
  WeightingFunction *m_wf;

  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Get the colour pair, energy integrators, etc.
   */
  virtual void setup();

public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  ThermostatLAWithWeightEC(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~ThermostatLAWithWeightEC() {}

  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase* p);
};

#endif
