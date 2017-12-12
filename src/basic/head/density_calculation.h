/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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
 * Add the particles self contribution to the local density
 * to the local density field.
 */
class DensityCalculation: public ParticleCache
{
 protected:
//   /*!
//   * Tag offset of the local density
   //   */
//  size_t m_offset;

   //DensityCalculation *pcl;
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
  * Initialise the property list
  */
  virtual void init();
  
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
   * LUT (2D Array), where all density values are stored
   */
  double **m_array_rho;

  /*!
   * Value of the interpolated local density.
   */
  double m_density_interpolation;

  /*!
   * Finds and approximate stored density values(LUT) with pressure and temperature values as input.
   * Bilinear interpolation is used for approximation.
   */
  virtual double calculateDensity(double inputT, double inputP, double m_Tmin, double m_pmin);

  /*!
   * Precalculates the density in given temperature and pressure ranges.
   * The density values are stored in a Look-Up table (2D Array) with fixed step sizes.
   */
  virtual void setupLUT(double m_Tmin, double m_pmin, double m_Tmax, double m_pmax,int arraysize_density, int arraysize_temperature);
  virtual void computeCacheFor(Particle* p) {
    double inputT = p->tag.doubleByOffset(m_temperatureOffset);
    double inputP = p->tag.doubleByOffset(m_pressureOffset);
    calculateDensity(inputT, inputP, m_Tmin, m_pmin);
    p->tag.doubleByOffset(m_offset) =  m_density_interpolation;
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
//     MSG_DEBUG("PressureCalculation::==", "called");
    if (typeid(c) == typeid(*this)) {

      /* DensityCalculation *cc = (DensityCalculation*) &c; */

      return true;
          /*m_wf->name() == cc->m_wf->name() && m_colour == cc->m_colour && m_stage == cc->m_stage && m_offset == cc->m_offset && m_symbolName == cc->m_symbolName;*/
    } else {
      return false;
    }
  } 
  /*!
   * If it belongs to a Node structure, setup this
   * ParticleCacheDensitySelfContribution
   */
  virtual void setup();


};

#endif
