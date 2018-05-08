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


#ifndef __PCA_IAPWS_IF97_H
#define __PCA_IAPWS_IF97_H

#ifdef HAVE_FREESTEAM

#include "particle_cache.h"

/*!
 * General \a ParticleCache computing a thermodynamic variable or a 
 * transport coefficient from two thermodynamic variables 'var1' and 
 * 'var2' representing a thermodynamic state. Calculations are 
 * performed based on IAPWS-IF97 (International Association for the 
 * Properties of Water and Steam. Revised release on the IAPWS 
 * industrial formulation 1997 for the thermodynamic properties of 
 * water and steam. adadad, August 2007) using the freesteam-library 
 * (http://freesteam.sourceforge.net/). Any instance of 
 * \a PCacheIAPWSIF97 or of a derived class uses a lookup table (LUT) 
 * in a predefined range 
 * [m_var1Min, m_var2Min] x [m_var1Max, m_var2Max] 
 * that is constructed once in advance to speedup the computation.
 */
class PCacheIAPWSIF97: public ParticleCache
{

 protected:

  /*!
   * The \a Symbol name of the first variable 'var1'
   */
  string m_var1Name;

  /*!
   * The \a Symbol name of the second variable 'var2'
   */
  string m_var2Name;

  /*!
   * \a Data tag offset of variable 'var1'
   */
  size_t m_var1Offset;

  /*!
   * \a Data tag offset of variable 'var2'
   */
  size_t m_var2Offset;

  /*!
   * Lower bound for variable 'var1' in the LUT 
   */
  double m_var1Min;

  /*!
   * Lower bound for variable 'var2' in the LUT 
   */
  double m_var2Min;

  /*!
   * Upper bound for variable 'var1' in the LUT 
   */
  double m_var1Max;

  /*!
   * Upper bound for variable 'var2' in the LUT 
   */
  double m_var2Max;

  /*!
   * Number of look-up values for variable 'var1' in the LUT.
   */
  size_t m_arraySizeVar1;

  /*!
   * Number of look-up values for variable 'var2' in the LUT.
   */
  size_t m_arraySizeVar2;

  /*!
   * Step size of look-up values for variable 'var2' internally derived 
   * from \a m_var1Min, \a m_var1Max, and \a m_arraySizeVar1
   */
  double m_var1StepSize;

  /*!
   * Step size of look-up values for variable 'var2' internally derived 
   * from \a m_var2Min, \a m_var2Max, and \a m_arraySizeVar2
   */
  double m_var2StepSize;

  /*!
   * Look-up table (LUT, 2D Array), where all look-up values are stored
   */
  double *m_LUT;

  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   * Note the explicit copy of m_LUT
   */
  virtual ParticleCache* copyMySelf() {

    // The RHS should be called for an instantiatable subclass
    ParticleCache* tmpPC = this->shallowCopy();

    // points to same instance
    PCacheIAPWSIF97* tmpPCIAPWS = (PCacheIAPWSIF97*) tmpPC;

    size_t totSize = (m_arraySizeVar1)*(m_arraySizeVar2);
    
    tmpPCIAPWS->m_LUT = new double[totSize];

    for (size_t i = 0; i < totSize; ++i) {
      tmpPCIAPWS->m_LUT[i] = m_LUT[i];
    }
    
    return tmpPC;
  }
  
  /*!
   * Adds the expressions used by this \a Symbol to the given list.
   * Yes, this function does not have any symbols to add here, but this 
   * is an abstract class, so better make this function pure virtual.
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols) = 0;
  
  /*!
   * Returns the strings of those \a Symbols that the given class 
   * depends on due to hard-coded reasons (not due to runtime compiled 
   * expressions).
   * Here we know that any instantiatable subclass will have the two 
   * symbols to add as done below.
   * @param usedSymbols List to add the strings to.
   */
  virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
  {
    usedSymbols.push_back(m_var1Name);
    usedSymbols.push_back(m_var2Name);
  }
    
  /*!
   * Instances of subclasses check if their individual constraints are 
   * fulfilled
   */
  virtual void checkConstraints() = 0;

  /*!
   * Instances of subclasses return a new shallow copy by 
   * copy-constructor of their own class 
   */
  virtual ParticleCache* shallowCopy() = 0;
  
 public:
  /*!
   * Constructor
   * @param colour The particle's color
   * @param offset Tag offset of the local density
   * @param wf The weighting function to use for the local density calculation
   */
  PCacheIAPWSIF97
    (size_t colour, size_t offset, string symbolName);
  
  /*!
   * Constructor
   */
  PCacheIAPWSIF97
    (Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~PCacheIAPWSIF97(); 
  
  /*!
   * Precalculates the output value in the 'var1' and 'var2' ranges 
   * specified by the corresponding private member variables
   * The output values are stored in a Look-Up table (LUT, 2D Array) 
   * with fixed step sizes given by \a m_stepSizeVar1 and 
   * \a m_stepSizeVar2.
   */
  virtual void setupLUT();

  /*!
   * Calculates and stores the output for the given particle
   * @param p The given particle
   */
  virtual void computeCacheFor(Particle* p) {
    
    Data& pTag = p->tag;
    double& result = pTag.doubleByOffset(m_offset);
    
    calculateResult
      (
       result,
       pTag.doubleByOffset(m_var1Offset),
       pTag.doubleByOffset(m_var2Offset)
       ); 
  }

  /*!
   * In instances of subclasses, this function performs the actual
   * call of the appropriate freesteam function  
   * @param result Memory address for storing the result
   * @param inputVar1 Value of thermodynamic input variable 'var1'
   * @param inputVar2 Value of thermodynamic input variable 'var2'
   */
  virtual void freesteamCalculationForState(double& result, const double& inputVar1, const double& inputVar2) const = 0;

  /*!
   * General logic for determining the result by bilinear interpolation 
   * based on the 4 support values of the LUT which surround 
   * (\a inputVar1, \a inputVar2) 
   * @param result Memory address for storing the result
   * @param inputVar1 Value of thermodynamic input variable 'var1'
   * @param inputVar2 Value of thermodynamic input variable 'var2'
   */
  virtual void calculateResult(double& result, const double& inputVar1, const double& inputVar2) const;
  
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
   * \a PCacheIAPWSIF97
   */
  virtual void setup();

  /*!
   * Checks existence of input symbols required by this \a ParticleCache
   * in a hard-coded fashion (i.e., not through runtime compiled 
   * expressions).
   * @param colour Particle colour to be checked
   */
  virtual void checkInputSymbolExistences(size_t colour);
  
  /*!
   * Returns the values stored in the LUT
   */ 
  virtual double* returnLUTvals() {
    return m_LUT;
  }

  ////////// friends ///////////////
  
  /*!
   * Class that is allowed to set protected members for unittesting
   */
  friend class PCacheIAPWSIF97TestGetter;
  
};

#endif    // HAVE_FREESTEAM

#endif
