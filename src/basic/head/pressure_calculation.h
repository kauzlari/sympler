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


#ifndef __PRESSURE_CALCULATION_H
#define __PRESSURE_CALCULATION_H

#include "pca_iapws-if97.h"

/*!
 * Local pressure at the particle computed from the local density
 * and local temperature based on IAPWS-IF97 (International
 * Association for the Properties of Water and Steam. Revised
 * release on the IAPWS industrial formulation 1997 for the
 * thermodynamic properties of water and steam. adadad, August
 * 2007).
 */
class PressureCalculation: public PCacheIAPWSIF97
{
  
 protected:
  
  /*!
   * Initialise the property list
   */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   */
  virtual ParticleCache* shallowCopy()
  {
    return new PressureCalculation(*this);
  }
    
  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
  {
    
  }
  
  /*!
   * \a PressureCalculation checks here if a requested temperature is 
   * within region 3 (623.15K < T < 863.15K). See IAPWS-IF97 for more 
   * information
   */
  virtual void checkConstraints();


 public:
  /*!
   * Constructor
   * @param colour The particle's color
   * @param offset Tag offset of the local density
   * @param wf The weighting function to use for the local density calculation
   */
  PressureCalculation
    (size_t colour, size_t offset, string symbolName);
  
  /*!
   * Constructor
   */
  PressureCalculation
    (Simulation* parent);  

  /*!
   * Performs the actual call of the appropriate freesteam function 
   * which computes pressure from density (\a m_var1) and temperature 
   * (\a m_var2)  
   * @param result Memory address for storing the resulting pressure
   * @param inputVar1 Value of thermodynamic input variable 'var1'
   * (density)
   * @param inputVar2 Value of thermodynamic input variable 'var2'
   * (temperature)
   */
  virtual void freesteamCalculationForState
    (double& result, const double& inputVar1, const double& inputVar2)
    const;

  /*!
   * Take steps necessary to register this calculator
   */
  virtual void registerWithParticle();
  
  /*!
   * If it belongs to a Node structure, setup this instance of
   * \a PressureCalculation
   */
  virtual void setup();
  
};

#endif
