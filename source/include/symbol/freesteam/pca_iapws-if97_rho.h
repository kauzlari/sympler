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


#ifndef __PCA_IAPWSIF97_RHO_H
#define __PCA_IAPWSIF97_RHO_H

#ifdef HAVE_FREESTEAM

#include "pca_iapws-if97_2var.h"

/*!
 * Local density at the particle computed from the local pressure
 * and local temperature based on IAPWS-IF97 (International
 * Association for the Properties of Water and Steam. Revised
 * release on the IAPWS industrial formulation 1997 for the
 * thermodynamic properties of water and steam. adadad, August
 * 2007).
 */
class PCacheIAPWSIF97rho: public PCacheIAPWSIF97TwoVar
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
    return new PCacheIAPWSIF97rho(*this);
  }
    
  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
  {
    
  }

  /*!
   * \a PCacheIAPWSIF97rho does not require any checks, hence function
   * is empty
   */
  virtual void checkConstraints();


 public:
  /*!
   * Constructor
   * @param colour The particle's color
   * @param offset Tag offset of the local density
   * @param wf The weighting function to use for the local density calculation
   */
  PCacheIAPWSIF97rho
    (size_t colour, size_t offset, string symbolName);
  
  /*!
   * Constructor
   */
  PCacheIAPWSIF97rho
    (Simulation* parent);  

  /*!
   * Performs the actual call of the appropriate freesteam function 
   * which computes density from pressure (\a m_var1) and temperature 
   * (\a m_var2)  
   * @param result Memory address for storing the resulting density
   * @param inputVar1 Value of thermodynamic input variable 'var1'
   * (pressure)
   * @param inputVar2 Value of thermodynamic input variable 'var2'
   * (temperature)
   */
  virtual void freesteamCalculationForState(double& result) const;

  /*!
   * Take steps necessary to register this calculator
   */
  virtual void registerWithParticle();
  
  /*!
   * If it belongs to a Node structure, setup this instance of
   * \a PCacheIAPWSIF97rho
   */
  virtual void setup();
  
};

#endif    // HAVE_FREESTEAM

#endif
