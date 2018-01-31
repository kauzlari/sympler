/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


#ifndef __DENSITY_CALCULATION_H
#define __DENSITY_CALCULATION_H

#include "particle_cache.h"


/*!
 * Local density at the particle computed from the local pressure
 * and local temperature based on IAPWS-IF97 (International
 * Association for the Properties of Water and Steam. Revised
 * release on the IAPWS industrial formulation 1997 for the
 * thermodynamic properties of water and steam. adadad, August
 * 2007).
 */
class DensityCalculation: public ParticleCache
{
 protected:

  /*!
   * From which symbol to take the local temperature
   */
  string m_temperatureName;

  /*!
   * From which symbol to take the local pressure
   */
  string m_pressureName;

  /*!
   * From which symbol to take the local maximum temperature
   */
  string m_temperatureMaxName;

  /*!
   * From which symbol to take the local maximum pressure
   */
  string m_pressureMaxName;

  /*!
   * From which symbol to take the local minimum temperature
   */
  string m_temperatureMinName;

  /*!
   * From which symbol to take the local minimum pressure
   */
  string m_pressureMinName;

  /*!
   * Tag offset of the local temperature
   */
  size_t m_temperatureOffset;

  /*!
   * Tag offset of the local pressure
   */
  size_t m_pressureOffset;

  /*!
   * Tag offset of the local maximum temperature
   */
  size_t m_temperatureMaxOffset;

  /*!
   * Tag offset of the local maximum pressure
   */
  size_t m_pressureMaxOffset;

  /*!
   * Tag offset of the local minimum temperature
   */
  size_t m_temperatureMinOffset;

  /*!
   * Tag offset of the local minimum pressure
   */
  size_t m_pressureMinOffset;

  /*!
   * Value of the local minimum temperature
   */
  double m_Tmin;

  /*!
   * Value of the local minimum pressure
   */
  double m_pmin;

  /*!
   * Value of the local maximum temperature
   */
  double m_Tmax;

  /*!
   * Value of the local maximum pressure
   */
  double m_pmax;

  /*!
   * Temperature step size for calculating new density values.
   * The step size depends on the size of the LUT and the temperature range.
   */
  double m_calcstepT;

  /*!
   * Pressure step size for calculating new density values.
   * The step size depends on the size of the LUT and the pressure range.
   */
  double m_calcstepP;

  /*!
   * Size of the array for temperature values.
   */
  int m_arraysize_temperature;

  /*!
   * Size of the array for pressure values.
   */
  int m_arraysize_pressure;

  /*!
   * LUT (2D Array), where all density values are stored
   */
  double **m_array_rho;
  
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
  {
    
  }
  
  /*!
   * Returns the strings of those \a Symbols that the given class depends on
   * due to hard-coded reasons (not due to runtime compiled expressions).
   * @param usedSymbols List to add the strings to.
   */
  virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
  {
    usedSymbols.push_back(m_temperatureName);
    usedSymbols.push_back(m_pressureName);
  }

  
 public:

  /*!
   * Constructor
   * @param colour The particle's color
   * @param offset Tag offset of the local density
   * @param wf The weighting function to use for the local density calculation
   */
   DensityCalculation
       (size_t colour, size_t offset, string symbolName);
  
   /*!
    * Constructor
    */
  DensityCalculation
    (Simulation* parent);

  /*!
   * Destructor
   */
  virtual ~DensityCalculation();

  /*!
   * Finds and approximate stored density values(LUT) with pressure and temperature values as input.
   * Bilinear interpolation is used for approximation.
   */
  virtual double calculateDensity(double inputT, double inputP);

  /*!
   * Precalculates the density in given temperature and pressure ranges.
   * The density values are stored in a Look-Up table (2D Array) with fixed step sizes.
   */
  virtual void setupLUT();

  /*!
   * Calculates the density for the given particle
   * @param p The given particle
   */
  virtual void computeCacheFor(Particle* p) {

    Data& pTag = p->tag;
    
    pTag.doubleByOffset(m_offset) =
      calculateDensity(
		       pTag.doubleByOffset(m_temperatureOffset),
		       pTag.doubleByOffset(m_pressureOffset)
		       );
  }

  /*!
   * Take steps necessary to register this calculator
   */
  virtual void registerWithParticle();

  /*!
   * Does this calculator equal \a c?
   * @param c Other calculator
   */
  virtual bool operator==(const ParticleCache &c) const {

    if (typeid(c) == typeid(*this)) {
      return true;
          /*m_wf->name() == cc->m_wf->name() && m_colour == cc->m_colour && m_stage == cc->m_stage && m_offset == cc->m_offset && m_symbolName == cc->m_symbolName;*/
    } else {
      return false;
    }
  }
  
  /*!
   * If it belongs to a Node structure, setup this instance of
   * \a DensityCalculation
   */
  virtual void setup();

  /*!
   * Returns the values stored in the LUT
   */ 
  virtual double** returnLUTvals() {
    return m_array_rho;
  }

};

#endif
