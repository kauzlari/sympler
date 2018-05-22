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

#include "fake_function_particle.h"

#include <cmath>
#include <limits>

CPPUNIT_TEST_SUITE_REGISTRATION (SymbolFParticleVelsTest);

void SymbolFParticleVelsTest :: setUp (void)
{

  // This is required since SymbolFParticleVels::setup() calls the
  // Controller of Simulation
  // This starts getting bulky, but the pain is still below
  // threshold...
  Simulation* sim = new Simulation();
  Controller* cont = new Controller(sim);
  sim -> setController(cont);
  
  // For the tests implemented so far (2018-04-20), NULL is fine
  m_symbolFParticleVels = new SymbolFParticleVels(/* Simulation* */ sim);
  m_symbolFParticleVels -> setFunctionParticle(new FakeFunctionParticle());
}

void SymbolFParticleVelsTest :: tearDown (void) 
{
  delete m_symbolFParticleVels;
}

void SymbolFParticleVelsTest :: initTest (void)
{
  // NOTE: SymbolFParticleVels::init() is called by the constructor
  CPPUNIT_ASSERT_EQUAL (true, m_symbolFParticleVels -> overwriting());
  CPPUNIT_ASSERT_EQUAL (string("__force__"), m_symbolFParticleVels -> mySymbolName());
}

void SymbolFParticleVelsTest :: setupTest (void)
{
  string exceptMsg = "ERROR: Unspecified exception in SymbolFParticleVelsTest :: setupTest for module " + m_symbolFParticleVels->className() + ". Please check!";

  // We catch the first gError and do nothing, since we do not yet
  // (2018-04-23) test xml-input. FIXME: could we do that with some
  // FakeClass?
  try {
    m_symbolFParticleVels -> setup();
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }

  // Hence, currently (2018-04-23) we expect exceptions
  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg, m_symbolFParticleVels -> setup(), gError);
  
  // Yes we repeat it, because some setups might have messed things up.
  // FIXME: Currently (2018-04-23) the called SymbolFParticleVels::setup()
  // does not do much, since it immediately stops after throwing an
  // exception. Probably this can only be fixed by adding tests for
  // the whole ParticleCache hierarchy
  SymbolFParticleVelsTest :: initTest ();
}

void SymbolFParticleVelsTest :: setupOffsetTest (void)
{
  try {
    m_symbolFParticleVels -> setup();
  } catch(gError& err) {}

  // FIXME: This is currently (2018-04-23) not a true test since
  // calling m_symbolFParticleVels -> setupOffset() is protected and setup()
  // throws an exception before arriving at setupOffset(). the
  // assertion will pass since we do m_offset = 2nd argument  in the
  // constructor of SymbolFParticleVels as well
  CPPUNIT_ASSERT_EQUAL (m_symbolFParticleVels -> offset(),
			std::numeric_limits<std::size_t>::max());
}

void SymbolFParticleVelsTest :: computeCacheForTest (void)
{
  
  Particle* p = new Particle();

  // first force slot used
  m_symbolFParticleVels -> setForceIndexTo(0);
  m_symbolFParticleVels -> computeCacheFor(p);
  // last force slot used
  m_symbolFParticleVels -> setForceIndexTo(FORCE_HIST_SIZE-1);
  m_symbolFParticleVels -> computeCacheFor(p);

  // looks a bit pointless to me that test...
  for(size_t i = 0; i < SPACE_DIMS; ++i) {
    // that's what FakeFunctionParticle::computeCacheFor(p) should have
    // done:
    CPPUNIT_ASSERT_EQUAL ((double) i, (p->force[0])[i]);
    CPPUNIT_ASSERT_EQUAL ((double) i, (p->force[FORCE_HIST_SIZE-1])[i]);

  }
}

