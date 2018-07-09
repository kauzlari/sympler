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



#ifndef __PCA_IAPWSIF97_1VAR_TEST_H
#define __PCA_IAPWSIF97_1VAR_TEST_H

#ifdef HAVE_FREESTEAM

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "pca_iapws-if97_1var.h"

using namespace std;


/*!
 * Class that is allowed to set protected members of class 
 * \a PCacheIAPWSIF97OneVar since \a PCacheIAPWSIF97OneVarTestGetter is a friend of 
 * \a PCacheIAPWSIF97OneVar.
 */
class PCacheIAPWSIF97OneVarTestGetter
{
  
 public: 

  /*!
   * Constructor
   * @param pc Instance of \a PCacheIAPWSIF97OneVar to be set up
   */
  PCacheIAPWSIF97OneVarTestGetter(PCacheIAPWSIF97OneVar* pc)
    : m_pc(pc)
    {}
  
  /*!
   * Destructor
   */
  virtual ~PCacheIAPWSIF97OneVarTestGetter() {
  }
  
  /*!
   * Set protected member \a PCacheIAPWSIF97OneVar::m_var1Max to given 
   * value
   * @param var1Max Value that m_var1Max should be set to
   */
  void set_m_var1Max(const double& var1Max) {
    m_pc->m_var1Max = var1Max;
  }
    
  /*!
   * Set protected member \a PCacheIAPWSIF97OneVar::m_var1Min to given 
   * value
   * @param var1Min Value that m_var1Min should be set to
   */
  void set_m_var1Min(const double& var1Min) {
    m_pc->m_var1Min = var1Min;
  }
    
  /*!
   * Set protected member \a PCacheIAPWSIF97OneVar::m_arraySizeVar1 to given 
   * value
   * @param arraysizeVar1 Value that m_arraySizeVar1 should be set to
   */
  void set_m_arraySizeVar1(const double& arraySizeVar1) {
    m_pc->m_arraySizeVar1 = arraySizeVar1;
  }

  /*!
   * Return protected member \a PCacheIAPWSIF97OneVar::m_inputVarPtrs 
   * Since it is a pointer, it also allows to set values in the 
   * original array
   * @param arraysizeVar1 Value that m_arraySizeVar1 should be set to
   */
  vector<double*>* return_m_inputVarPtrs() {
    return &(m_pc->m_inputVarPtrs);
  }

  /*!
   * Call \a PCacheIAPWSIF97OneVar::copyMySelf() or of childclass
   */
  ParticleCache* get_m_pc_copyMySelf() {
    return m_pc->copyMySelf();
  }
    

 protected:
  
  PCacheIAPWSIF97OneVar* m_pc;
  
}; // end of PCacheIAPWSIF97OneVarTestGetter


/*!
 * Abstract parent test class for abstract parent class 
 * \a PCacheIAPWSIF97OneVar
 */
class PCacheIAPWSIF97OneVarTest : public CPPUNIT_NS :: TestFixture
{
  
 public:

  virtual void setUp (void);
  virtual void tearDown (void);

 protected:

  virtual void setupLUTTest (void);

  /*!
   * Tests PCacheIAPWSIF97OneVar::calculateResult(..)
   * Implemented by instantiatable subclasses of PCacheIAPWSIF97OneVarTest
   * which are testing instantiatable subclasses of PCacheIAPWSIF97OneVar
   */
  virtual void calculateResultTest (void) = 0;

  virtual void copyMySelfTest();
  
  Simulation *m_simulation; 
  PCacheIAPWSIF97OneVar *m_pc;
  PCacheIAPWSIF97OneVarTestGetter *m_pcGetter;
  
  double m_var1Max;
  double m_var1Min;
  size_t m_arraySizeVar1;

  /*!
   * Helper for multiple execution
   * @param var1 First variable describing the thermodynamic state
   */
  virtual void execCalculationTest(const double& var1);
  
};

#endif    // HAVE_FREESTEAM

#endif
