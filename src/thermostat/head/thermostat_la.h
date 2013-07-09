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




#ifndef __THERMOSTAT_LA_H
#define __THERMOSTAT_LA_H 


#include "random.h"
#include "simulation.h"
#include "general.h"
#include "callable.h"
#include "colour_pair.h"
#include "thermostat_with_rng.h"
#include "function_pair.h"


using namespace std;

class Phase;
class Simulation;


class ThermostatLA : public ThermostatWithRng
{
protected:
  /*!
   * The pair of colors to thermalize
   */
  ColourPair* m_cp;

  /*!
   * The species to thermalzie
   */
  pair<string, string> m_species;

  /*!
   * The cut-off radius
   */
  double m_cutoff;

  /*!
   * Thermalization probability = dissipation * time step
   */
  FunctionPair m_probability;
/*   double m_probability; */

  /*!
   * First particle factor multiplying the relative velocity drawn by the thermostat
   */
  FunctionPair m_firstFac;

  /*!
   * Second particle factor multiplying the relative velocity drawn by the thermostat
   */
  FunctionPair m_secondFac;

  /*!
   * First particle addend added to the relative velocity drawn by the thermostat
   */
  FunctionPair m_firstAdd;

  /*!
   * Second particle addend added to the relative velocity drawn by the thermostat
   */
  FunctionPair m_secondAdd;

  /*!
   * The temperature of the heat bath
   */
  double m_temperature;

  /*!
   * Width of the equilibrium distribution
   */
  double m_MaxwellFactor;
			
  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  ThermostatLA(Simulation *sim);

  /*!
   * Destructor
   */
  virtual ~ThermostatLA();

  /*!
   * Get the \a ColourPair and initialize variables
   */
  virtual void setup();

  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase *p);
};


#endif

