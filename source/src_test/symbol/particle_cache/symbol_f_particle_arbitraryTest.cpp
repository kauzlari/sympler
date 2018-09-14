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


#include "symbol_f_particle_arbitraryTest.h"

#include "fake_function_particle.h"

#include <cmath>
#include <limits>

// commented out since this is an abstract class testing an abstract
// class through inheriting children
// CPPUNIT_TEST_SUITE_REGISTRATION (SymbolFParticleArbitraryTest);

void SymbolFParticleArbitraryTest :: setUp (void)
{
  // This is required since SymbolFParticleArbitrary::setup() calls the
  // Controller of Simulation. This starts getting bulky, but the pain
  // is still below threshold...
  m_simulation = new Simulation();
  Controller* cont = new Controller(m_simulation);
  m_simulation -> setController(cont);  
}

// helper function; not a test; meant to be called by setUp of child test class
void SymbolFParticleArbitraryTest :: setUpFunctionParticle() {
  
  m_symbolFParticleArbitrary -> setFunctionParticle(new FakeFunctionParticle());
}


void SymbolFParticleArbitraryTest :: tearDown (void) 
{
	delete m_simulation -> controller();
  delete m_simulation;
	delete m_symbolFParticleArbitrary;
}

void SymbolFParticleArbitraryTest :: initTest (void)
{
  // NOTE: SymbolFParticleArbitrary::init() is called by the constructor
  CPPUNIT_ASSERT_EQUAL
    (true, m_symbolFParticleArbitrary -> overwriting());
  CPPUNIT_ASSERT_EQUAL
    (string("PREDEFINED"), m_symbolFParticleArbitrary -> mySymbolName());
}

void SymbolFParticleArbitraryTest :: setupTest (void)
{

  string exceptMsg = "ERROR: Unspecified exception in SymbolFParticleArbitraryTest :: setupTest for module " + m_symbolFParticleArbitrary->className() + ". Please check!";

  // We catch the first gError and do nothing else, since we do not yet
  // (2018-06-07) test xml-input.
  // FIXME: could we do that with some FakeClass?
  // FIXME: registration for precomputation not tested yet
  try {
    // better use a copy such that setup() is not run twice
    SymbolFParticleArbitrary* copy =
      (SymbolFParticleArbitrary*) (m_symbolFParticleArbitrary -> returnCopy());
    copy -> setup();
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }

  // Hence, currently (2018-06-07) we expect exceptions (thrown by
  // ParticleCacheArbitrary or grand parents) due to unset xml-input
  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg, m_symbolFParticleArbitrary -> setup(), gError);

}

void SymbolFParticleArbitraryTest :: setupOffsetTest (void)
{
  // exceptions will be thrown but we are not interested in them in
  // this test...
  try {
    m_symbolFParticleArbitrary -> setup();
  } catch(gError& err) {/* ...so we do nothing */}

  // FIXME: the right thing is tested, but this is currently
  // (2018-06-07) not a true test since
  // calling m_symbolFParticleArbitrary -> setupOffset() is protected
  // and setup() throws an exception before arriving at setupOffset().
  // The assertion will pass since we do m_offset = 2nd argument below
  // in the constructor of SymbolFParticleArbitrary as well
  CPPUNIT_ASSERT_EQUAL (m_symbolFParticleArbitrary -> offset(),
			std::numeric_limits<std::size_t>::max());
}

void SymbolFParticleArbitraryTest :: precomputeTest (void)
{
  m_symbolFParticleArbitrary -> precompute();

  CPPUNIT_ASSERT_EQUAL
    (m_symbolFParticleArbitrary -> offset(),
     (size_t) m_simulation -> controller() -> forceIndex()
     );
}

void SymbolFParticleArbitraryTest :: setForceIndexToTest (void)
{
  size_t testNum = 12345;
  
  m_symbolFParticleArbitrary -> setForceIndexTo(testNum);

  CPPUNIT_ASSERT_EQUAL
    (m_symbolFParticleArbitrary -> offset(), testNum);
}

void SymbolFParticleArbitraryTest :: computeCacheFor(Particle* p, size_t forceIndex) {

  m_symbolFParticleArbitrary -> setForceIndexTo(forceIndex);

  m_symbolFParticleArbitrary -> computeCacheFor(p);

}
