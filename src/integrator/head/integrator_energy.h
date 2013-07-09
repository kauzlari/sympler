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



#ifndef __INTEGRATOR_ENERGY_H
#define __INTEGRATOR_ENERGY_H 

#include <math.h>

#include "function_fixed.h"
#include "particle.h"
#include "integrator_scalar.h"

using namespace std;

class GenF;
class Phase;
class Controller;

/*!
 * Integrator for the energy degree of freedom that can calculate the temperature
 * for each particle
 */
class IntegratorEnergy: public IntegratorScalar
{
protected:
  /*!
   * Function giving the temperature for each particle
   */
  FunctionFixed m_ds_de;

  /*!
   * Tag offset of the temperature
   */
  size_t m_T_offset;

  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Return the reciprocal of the temperature
   * @param energy Calculate the reciprocal temperature for this energy
   */
  double dEntropy_dEnergy(double energy) {
    return m_ds_de(energy);
  }

public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorEnergy(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorEnergy();

  /*!
   * Register the field and the force of the field with the \a Particle
   */
  virtual void setup();

  /*!
   * Calculate and store the temperature for each particle
   */
  virtual void deriveQuantities();

  /*!
   * Return the tag offset for the temperature
   */
  size_t TOffset() const {
    return m_T_offset;
  }

#ifdef _OPENMP
  virtual string dofIntegr();


#endif

  /*!
   * Return the reciprocal of the temperature
   * @param energy Calculate the temperature for this energy
   */
  double reciprocalTemperature(double energy) {
    return dEntropy_dEnergy(energy);
  }

  /*!
   * Return the reciprocal of the temperature
   * @param p Calculate the temperature for this \a Particle
   */
  double reciprocalTemperature(Particle &p) {
    return dEntropy_dEnergy
      (p.tag.doubleByOffset(m_scalar_offset)
       );
  }

  /*!
   * Return the temperature
   * @param energy Calculate the temperature for this energy
   */
  double temperature(double energy) {
    return 1/dEntropy_dEnergy(energy);
  }

  /*!
   * Return the temperature
   * @param p Calculate the temperature for this \a Particle
   */
  double temperature(Particle &p) {
    return 1/reciprocalTemperature(p);
  }
};

#endif
