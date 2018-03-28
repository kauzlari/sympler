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


#include "pca_iapws-if97Test.h"

#include "simulation.h"


void PCacheIAPWSIF97Test  :: setUp (void)
{

  // Size of the arrays
  m_arraySizeVar1 = 10;
  m_arraySizeVar2 = 10;
  
  // Arguments must have been set in setUp of subclass before
  m_pcGetter = new PCacheIAPWSIF97TestGetter(m_pc);

  m_pcGetter->set_m_var1Max(m_var1Max);
  m_pcGetter->set_m_var2Max(m_var2Max);
  m_pcGetter->set_m_var1Min(m_var1Min);
  m_pcGetter->set_m_var2Min(m_var2Min);
  m_pcGetter->set_m_arraySizeVar1(m_arraySizeVar1);
  m_pcGetter->set_m_arraySizeVar2(m_arraySizeVar2);
  
  m_pc->setupLUT();

}


/*!
 * Delete objects
 */
void PCacheIAPWSIF97Test  :: tearDown (void) 
{
  delete m_simulation; 
  delete m_pcGetter;
  delete m_pc;
}

void PCacheIAPWSIF97Test  :: setupLUTTest (void) {

  // freesteamCalculationForState(..) was already called in setUp(..)
  // when calling m_pc->setupLUT(..), so exceptions were already thrown
  // if necessary. Hence, we only have to check the LUT here for some
  // correct entries.

  double result;

  // test of minimum
  m_pc -> freesteamCalculationForState(result, m_var1Min, m_var2Min);
  CPPUNIT_ASSERT_EQUAL (m_pc -> returnLUTvals()[0], result);

  // test of maximum
  m_pc -> freesteamCalculationForState(result, m_var1Max, m_var2Max);
  CPPUNIT_ASSERT_EQUAL (m_pc -> returnLUTvals()[m_arraySizeVar1*m_arraySizeVar2 - 1], result);

}


void PCacheIAPWSIF97Test  :: copyMySelfTest (void) {

  PCacheIAPWSIF97* tmpPCIAPWS = (PCacheIAPWSIF97*)(m_pcGetter->get_m_pc_copyMySelf());

  // double* pointer should NOT be equal
  CPPUNIT_ASSERT_EQUAL (tmpPCIAPWS->returnLUTvals() != m_pc->returnLUTvals(), true);
  for(size_t i = 0; i < m_arraySizeVar1; ++i) {
    size_t slot = i*m_arraySizeVar2;
    for(size_t j = 0; j < m_arraySizeVar2; ++j) {
      // values should be equal
      CPPUNIT_ASSERT_EQUAL (tmpPCIAPWS->returnLUTvals()[slot], m_pc->returnLUTvals()[slot]);
      ++slot;
    }
  } 
}

void PCacheIAPWSIF97Test  :: execCalculationTest(const double& var1, const double& var2) {
  
  // Assertion to check if the interpolation is correct.
  double interpol;
  double result;
  string exceptMsg;
  
  exceptMsg = "ERROR: Unspecified exception in PCacheIAPWSIF97Test::execCalculationTest for module " + m_pc->className() + ". Please check!";
  m_pc -> freesteamCalculationForState(result, var1, var2);
  try {
    m_pc -> calculateResult(interpol, var1, var2);
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }
  CPPUNIT_ASSERT_NO_THROW_MESSAGE(exceptMsg, m_pc -> calculateResult(interpol, var1, var2));
  CPPUNIT_ASSERT_DOUBLES_EQUAL (interpol, result, interpol/100.); // relativ delta: 1% of expected
  
}
