/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#include "symbol_f_particle_scalarTest.h"

#include "fake_symbol_f_particle_scalar.h"

#include <cmath>
#include <limits>

CPPUNIT_TEST_SUITE_REGISTRATION (SymbolFParticleScalarTest);


void SymbolFParticleScalarTest :: setUp (void)
{
  SymbolFParticleArbitraryTest::setUp();
  
  // m_simulation is assumed to be set in
  // SymbolFParticleArbitraryTest::setUp()
  // Search for "FakeSymbolFParticleScalar" below to find cases where we
  // really use the FakeClass properties
  m_symbolFParticleArbitrary = new FakeSymbolFParticleScalar(m_simulation);

  ((FakeSymbolFParticleScalar*) m_symbolFParticleArbitrary)
    -> setColour(0);
  
  setUpFunctionParticle();
}


void SymbolFParticleScalarTest :: initTest (void)
{
  SymbolFParticleArbitraryTest :: initTest();

  CPPUNIT_ASSERT_EQUAL
    (string("UNDEFINED"),
     ((SymbolFParticleScalar*) m_symbolFParticleArbitrary)
     -> returnVariableName()
     );  
}


void SymbolFParticleScalarTest :: setupTest (void)
{

  // here, the actually to be tested setup() will be called
  SymbolFParticleArbitraryTest :: setupTest();

  string variableName =
    ((SymbolFParticleScalar*) m_symbolFParticleArbitrary)
    -> returnVariableName();
  
  CPPUNIT_ASSERT_EQUAL
    (string("force_" + variableName),
     m_symbolFParticleArbitrary -> mySymbolName()
     );
}


void SymbolFParticleScalarTest :: computeCacheForTest (void)
{
  
  string variableName =
    ((SymbolFParticleScalar*) m_symbolFParticleArbitrary)
    -> returnVariableName();

  string forceName = string("TESTforce_" + variableName + "_");
  
  
  size_t colour = 0;
  
  Particle::s_tag_format.resize(1);
  Particle::s_tag_format[colour].setFormatAndAlloc(new DataFormat());

  // FIXME: ideally we would get the force offsets directly from SymbolFParticleScalar as a further test. But the required method setupOffset() cannot yet be run within these unittests due to thrown exceptions. Hence we set the offsets manually to the return value of the manually called Particle::s_tag_format[colour].addAttribute(..)
  // size_t firstForceOffset =
  //   ((SymbolFParticleScalar*) m_symbolFParticleArbitrary)
  //   -> returnForceOffset(0);

  // size_t lastForceOffset =
  //   ((SymbolFParticleScalar*) m_symbolFParticleArbitrary)
  //   -> returnForceOffset(FORCE_HIST_SIZE-1);

  // this reallocs memory if necessary
  size_t firstForceOffset
    = Particle::s_tag_format[colour].addAttribute
    (forceName + "0", DataFormat::DOUBLE, false, forceName + "0").offset;

  // last force offset for general FORCE_HIST_SIZE
  size_t lastForceOffset;
  for (size_t i = 1; i < FORCE_HIST_SIZE; ++i) {
    lastForceOffset = Particle::s_tag_format[colour].addAttribute(forceName + ObjToString(FORCE_HIST_SIZE-1), DataFormat::DOUBLE, false, forceName + ObjToString(FORCE_HIST_SIZE-1)).offset;
  }

  // with the chosen names above, SymbolFParticleScalar should find the
  // variables in the tag and the same offsets as in firstForceOffset
  // and lastForceOffset
  ((FakeSymbolFParticleScalar*) m_symbolFParticleArbitrary)
    -> setupForceOffsetForTestsOnly();
    
  
  // creation  for specific colour ensures reservation of memory for
  // tag as defined by Particle::s_tag_format[colour]
  Particle* p = new Particle(colour);
  
  // 2nd argument is force index. The called function calls the actually
  // to be tested function
  // first force slot used
  SymbolFParticleArbitraryTest :: computeCacheFor(p, 0 /*force index*/);

  // last force slot used
  SymbolFParticleArbitraryTest
    :: computeCacheFor(p, FORCE_HIST_SIZE-1 /*force index*/);
  
  // looks a bit pointless to me that test...
  // ... but that's what FakeFunctionParticle::computeCacheFor(p)
  // should have done according to our own fake-definition:  
  CPPUNIT_ASSERT_EQUAL
    ( 3., p -> tag.doubleByOffset(firstForceOffset) );
  
  CPPUNIT_ASSERT_EQUAL
    ( 3., p -> tag.doubleByOffset(lastForceOffset) );    

  delete p;  

}
