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



#ifndef __PCA_IAPWSIF97_2VAR_TEST_H
#define __PCA_IAPWSIF97_2VAR_TEST_H

#ifdef HAVE_FREESTEAM

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "pca_iapws-if97_2var.h"
#include "pca_iapws-if97_1varTest.h"

using namespace std;


/*!
 * Class that is allowed to set protected members of class 
 * \a PCacheIAPWSIF97 since \a PCacheIAPWSIF97TestGetter is a friend of 
 * \a PCacheIAPWSIF97.
 */
class PCacheIAPWSIF97TwoVarTestGetter : public PCacheIAPWSIF97OneVarTestGetter
{
  
 public: 

  /*!
   * Constructor
   * @param pc Instance of \a PCacheIAPWSIF97OneVar or subclass to be 
   * set up
   */
  PCacheIAPWSIF97TwoVarTestGetter(PCacheIAPWSIF97OneVar* pc)
    : PCacheIAPWSIF97OneVarTestGetter(pc)
    {}
  
  /*!
   * Destructor
   */
  virtual ~PCacheIAPWSIF97TwoVarTestGetter() {
  }
  
  /*!
   * Set protected member \a PCacheIAPWSIF97TwoVar::m_var2Max to given 
   * value
   * @param var2Max Value that m_var2Max should be set to
   */
  void set_m_var2Max(const double& var2Max) {
    // type cast to access PCacheIAPWSIF97TwoVar-specific method
    ((PCacheIAPWSIF97TwoVar*) m_pc)->m_var2Max = var2Max;
  }
  
  /*!
   * Set protected member \a PCacheIAPWSIF97TwoVar::m_var2Min to given 
   * value
   * @param var2Min Value that m_var2Min should be set to
   */
  void set_m_var2Min(const double& var2Min) {
    // type cast to access PCacheIAPWSIF97TwoVar-specific method
    ((PCacheIAPWSIF97TwoVar*) m_pc)->m_var2Min = var2Min;
  }
  
  /*!
   * Set protected member \a PCacheIAPWSIF97TwoVar::m_arraySizeVar2 to given 
   * value
   * @param arraysizeVar2 Value that m_arraySizeVar2 should be set to
   */
  void set_m_arraySizeVar2(const double& arraySizeVar2) {
    // type cast to access PCacheIAPWSIF97TwoVar-specific method
    ((PCacheIAPWSIF97TwoVar*) m_pc)->m_arraySizeVar2 = arraySizeVar2;
  }
    

 protected:
  

}; // end of PCacheIAPWSIF97TwoVarTestGetter


/*!
 * Abstract parent test class for abstract parent class 
 * \a PCacheIAPWSIF97TwoVar
 */
class PCacheIAPWSIF97TwoVarTest : public PCacheIAPWSIF97OneVarTest
{
  
 public:

  /*!
   * Note: For simplicity this function currently completely overwrites 
   * the parent class function
   */
  virtual void setUp (void);

  /* virtual void tearDown (void); */

 protected:

  virtual void setupLUTTest (void);

  /*!
   * Tests PCacheIAPWSIF97TwoVar::calculateResult(..)
   * Implemented by instantiatable subclasses of PCacheIAPWSIF97TwoVarTest
   * which are testing instantiatable subclasses of PCacheIAPWSIF97TwoVar
   */
  virtual void calculateResultTest (void) = 0;

  virtual void copyMySelfTest();
  
  /* Simulation *m_simulation;  */
  /* PCacheIAPWSIF97TwoVar *m_pc; */
  /* PCacheIAPWSIF97TwoVarTestGetter *m_pcGetter; */
  
  double m_var2Max;
  double m_var2Min;
  size_t m_arraySizeVar2;

  /*!
   * Helper for multiple execution, inherited from parent class. For
   * now (2018-04-04), we let it for simplicity throw an exception, 
   * since we must use two variables here, but we do not want to 
   * introduce again a vector-member as in the tested class. 
   * @param var1 First variable describing the thermodynamic state
   */
  virtual void execCalculationTest(const double& var1) {
    throw gError("PCacheIAPWSIF97TwoVarTest::execCalculationTest(var1)"
		 , "should not have been called for a child class of "
		 "PCacheIAPWSIF97TwoVarTest. Something went wrong in "
		 "the implementation of a test class. Aborting.");
  }
  
  /*!
   * Helper for multiple execution. This method is newly introduced in
   * this class and is not inherited from parent classes. For now 
   * (2018-04-04) this is a simple soluton for the fact that we must 
   * use two variables here.
   * @param var1 First variable describing the thermodynamic state
   * @param var2 Second variable describing the thermodynamic state
   */
  virtual void execCalculationTest(const double& var1, const double& var2);
  
};

#endif    // HAVE_FREESTEAM

#endif
