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


#include "symbol_f_particle_velsTest.h"

#include <cmath>
#include <limits>

CPPUNIT_TEST_SUITE_REGISTRATION (SymbolFParticleVelsTest);


void SymbolFParticleVelsTest :: setUp (void)
{
  SymbolFParticleArbitraryTest::setUp();
  
  // m_simulation assumed to be set in SymbolFParticleArbitraryTest::setUp()
  m_symbolFParticleArbitrary = new SymbolFParticleVels(m_simulation);

  setUpFunctionParticle();
}


void SymbolFParticleVelsTest :: setupTest (void)
{
  // here, the actually to be tested setup() will be called
  SymbolFParticleArbitraryTest :: setupTest();
    
  CPPUNIT_ASSERT_EQUAL (string("__force__"), m_symbolFParticleArbitrary -> mySymbolName());
}


void SymbolFParticleVelsTest :: computeCacheForTest (void)
{  
  Particle* p = new Particle();

  // 2nd argument is force index. The called function calls the actually
  // to be tested function
  // first force slot used
  SymbolFParticleArbitraryTest :: computeCacheFor(p, 0);
  // last force slot used
  SymbolFParticleArbitraryTest :: computeCacheFor(p, FORCE_HIST_SIZE-1);
  
  // looks a bit pointless to me that test...
  for(size_t i = 0; i < SPACE_DIMS; ++i) {
    // that's what FakeFunctionParticle::computeCacheFor(p) should have
    // done according to our own fake-definition:
    CPPUNIT_ASSERT_EQUAL ((double) i, (p->force[0])[i]);
    CPPUNIT_ASSERT_EQUAL ((double) i, (p->force[FORCE_HIST_SIZE-1])[i]);    
  }

  delete p;  
}

