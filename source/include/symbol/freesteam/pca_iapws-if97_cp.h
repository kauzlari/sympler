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


#ifndef __PCA_IAPWSIF97_CP_H
#define __PCA_IAPWSIF97_CP_H

#ifdef HAVE_FREESTEAM

#include "pca_iapws-if97_1var.h"

/*!
 * Local isobaric heat capacity at the particle computed from the local 
 * temperature at a user-defined constant pressure
 * based on IAPWS-IF97 (International Association for the Properties of 
 * Water and Steam. Revised release on the IAPWS industrial formulation 
 * 1997 for the thermodynamic properties of water and steam. adadad, 
 * August 2007).
 */
class PCacheIAPWSIF97Cp: public PCacheIAPWSIF97OneVar
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
    return new PCacheIAPWSIF97Cp(*this);
  }
    
  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
  {
    
  }

  /*!
   * \a PCacheIAPWSIF97Cp does not require any checks, hence function
   * is empty
   */
  virtual void checkConstraints();

  /*!
   * Value of the constant parameter
   * FIXME: Here, for the isobaric heat capacity (Cp), it is pressure. 
   * Check if this is worth a parent class when, e.g., Cv is added.
   * Probably not.
   */
  double m_constP;
  
 public:
  /*!
   * Constructor
   * @param colour The particle's color
   * @param offset Tag offset of the local density
   * @param wf The weighting function to use for the local density calculation
   */
  PCacheIAPWSIF97Cp
    (size_t colour, size_t offset, string symbolName);
  
  /*!
   * Constructor
   */
  PCacheIAPWSIF97Cp
    (Simulation* parent);  

  /*!
   * Performs the actual call of the appropriate freesteam function 
   * which computes isobaric heat capacity from temperature 
   * (\a m_inputVarPtrs) at a given constant pressure (\a m_constP). 
   * Note the mentioned class members in brackets which make input 
   * arguments obsolete.  
   * @param result Memory address for storing the resulting density
   */
  /* docu of obsolete arguments
   * @param inputVar1 Value of thermodynamic input variable 'var1'
   * (density)
   * @param inputVar2 Value of thermodynamic input variable 'var2'
   * (temperature)
   */
  virtual void freesteamCalculationForState
    (double& result/* , const double& inputVar1, const double& inputVar2 */)
    const;

  /*!
   * Take steps necessary to register this calculator
   */
  virtual void registerWithParticle();
  
  /*!
   * If it belongs to a Node structure, setup this instance of
   * \a PCacheIAPWSIF97Cp
   */
  virtual void setup();
  
  /*!
   * 3 characters defining the derivative dz/dx|y=const with 
   * thermodynamic variables to be taken. This array is required as 
   * input by the freesteam function. 
   * freesteam_deriv(StreamState S, char[3] xyz).
   * For the computation of Cp, this is a never changing static and 
   * constant object defined in the cpp file.
   */
  static const char s_xyz[3];

  ////////// friends ///////////////
  
  /*!
   * Class that is allowed to set protected members for unittesting
   */
  friend class PCacheIAPWSIF97WithConstTestGetter;

};

#endif    // HAVE_FREESTEAM

#endif
