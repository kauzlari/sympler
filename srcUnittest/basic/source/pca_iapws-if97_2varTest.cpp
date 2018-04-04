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


#include "pca_iapws-if97_2varTest.h"

#include "simulation.h"


void PCacheIAPWSIF97TwoVarTest  :: setUp (void)
{

  // Size of the arrays
  m_arraySizeVar1 = 10;
  m_arraySizeVar2 = 10;
  
  // Arguments must have been set in setUp of subclass before
  m_pcGetter = new PCacheIAPWSIF97TwoVarTestGetter(m_pc);

  // handle to access the PCacheIAPWSIF97TwoVarTestGetter-specific
  // methods
  PCacheIAPWSIF97TwoVarTestGetter* m_pc2VarGetter
    = (PCacheIAPWSIF97TwoVarTestGetter*) m_pcGetter;
  
  m_pcGetter->set_m_var1Max(m_var1Max);
  m_pc2VarGetter->set_m_var2Max(m_var2Max);
  m_pcGetter->set_m_var1Min(m_var1Min);
  m_pc2VarGetter->set_m_var2Min(m_var2Min);
  m_pcGetter->set_m_arraySizeVar1(m_arraySizeVar1);
  m_pc2VarGetter->set_m_arraySizeVar2(m_arraySizeVar2);
  
  m_pc->setupLUT();

}


void PCacheIAPWSIF97TwoVarTest  :: setupLUTTest (void) {

  // freesteamCalculationForState(..) was already called in setUp(..)
  // when calling m_pc->setupLUT(..), so exceptions were already thrown
  // if necessary. Hence, we only have to check the LUT here for some
  // correct entries.

  // was vector<double*> PCacheIAPWSIF97OneVar::m_inputVarPtrs
  // correctly resized?
  // FIXME: this should in fact be part of a test of
  // PCacheIAPWSIF97OneVar::setup()
  CPPUNIT_ASSERT_EQUAL((*(m_pcGetter->return_m_inputVarPtrs())).size(), (size_t)2);
  
  // handle to double* PCacheIAPWSIF97OneVar::m_inputVarPtrs[0] for
  // setting the current input variable 1
  double** inputVar1 = &((*(m_pcGetter->return_m_inputVarPtrs()))[0]);

  // handle to double* PCacheIAPWSIF97OneVar::m_inputVarPtrs[1] for
  // setting the current input variable 2
  double** inputVar2 = &((*(m_pcGetter->return_m_inputVarPtrs()))[1]);

  double result;

  // test of minimum
  *inputVar1 = &m_var1Min;
  *inputVar2 = &m_var2Min;
  m_pc -> freesteamCalculationForState(result);
  CPPUNIT_ASSERT_EQUAL (m_pc -> returnLUTvals()[0], result);

  // test of maximum
  *inputVar1 = &m_var1Max;
  *inputVar2 = &m_var2Max;
  m_pc -> freesteamCalculationForState(result);
  // due to summation we get a slight error here, so far below 1/1e14
  CPPUNIT_ASSERT_DOUBLES_EQUAL (m_pc -> returnLUTvals()[m_arraySizeVar1*m_arraySizeVar2 - 1], result, result/1e14);

}


void PCacheIAPWSIF97TwoVarTest  :: copyMySelfTest (void) {

  PCacheIAPWSIF97TwoVar* tmpPCIAPWS = (PCacheIAPWSIF97TwoVar*)(m_pcGetter->get_m_pc_copyMySelf());

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

void PCacheIAPWSIF97TwoVarTest  :: execCalculationTest(const double& var1, const double& var2) {
  
  // Assertion to check if the interpolation is correct.
  double interpol;
  double result;
  string exceptMsg;

  // handle to double* PCacheIAPWSIF97OneVar::m_inputVarPtrs[0] for
  // setting the current input variable 1
  double** inputVar1 = &((*(m_pcGetter->return_m_inputVarPtrs()))[0]);
  *inputVar1 = (double*) &var1;

  // handle to double* PCacheIAPWSIF97OneVar::m_inputVarPtrs[1] for
  // setting the current input variable 2
  double** inputVar2 = &((*(m_pcGetter->return_m_inputVarPtrs()))[1]);
  *inputVar2 = (double*) &var2;
  
  exceptMsg = "ERROR: Unspecified exception in PCacheIAPWSIF97TwoVarTest::execCalculationTest for module " + m_pc->className() + ". Please check!";
  m_pc -> freesteamCalculationForState(result);
  try {
    m_pc -> calculateResult(interpol);
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }
  CPPUNIT_ASSERT_NO_THROW_MESSAGE(exceptMsg, m_pc -> calculateResult(interpol));
  CPPUNIT_ASSERT_DOUBLES_EQUAL (interpol, result, interpol/100.); // relativ delta: 1% of expected
  
}
