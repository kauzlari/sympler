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


#include "particle_velsTest.h"

#include "fake_function_particle.h"

#include <cmath>

CPPUNIT_TEST_SUITE_REGISTRATION (ParticleVelsTest);

void ParticleVelsTest :: setUp (void)
{
  // For the tests implemented so far (2018-04-20), NULL is fine
  m_particleVels = new ParticleVels(/* Simulation* or Node* */ NULL);
  m_particleVels -> setFunctionParticle(new FakeFunctionParticle());
}

void ParticleVelsTest :: tearDown (void) 
{
  delete m_particleVels;
}

void ParticleVelsTest :: initTest (void)
{
  // NOTE: ParticleVels::init() is called by the constructor
  CPPUNIT_ASSERT_EQUAL (m_particleVels -> overwriting(), true);
  CPPUNIT_ASSERT_EQUAL (m_particleVels -> mySymbolName(), string("v"));
}

void ParticleVelsTest :: setupTest (void)
{
  string exceptMsg = "ERROR: Unspecified exception in ParticleVelsTest :: setupTest for module " + m_particleVels->className() + ". Please check!";

  // We catch the first gError and do nothing, since we do not yet
  // (2018-04-23) test xml-input. FIXME: could we do that with some
  // FakeClass?
  try {
    m_particleVels -> setup();
  } catch(gError& err) {
    exceptMsg = err.message(); 
  }

  // Hence, currently (2018-04-23) we expect exceptions
  CPPUNIT_ASSERT_THROW_MESSAGE(exceptMsg, m_particleVels -> setup(), gError);
  
  // Yes we repeat it, because some setups might have messed things up.
  // FIXME: Currently (2018-04-23) the called ParticleVels::setup()
  // does not do much, since it immediately stops after throwing an
  // exception. Probably this can only be fixed by adding tests for
  // the whole ParticleCache hierarchy
  ParticleVelsTest :: initTest ();
}

void ParticleVelsTest :: setupOffsetTest (void)
{
  try {
    m_particleVels -> setup();
  } catch(gError& err) {}

  // FIXME: This is currently (2018-04-23) not a true test since
  // calling m_particleVels -> setupOffset() is protected and setup()
  // throws an exception before arriving at setupOffset(). the
  // assertion will pass since we do m_offset = HUGE_VAL  in the
  // constructor of ParticleVels as well
  CPPUNIT_ASSERT_EQUAL (m_particleVels -> offset(), (size_t) HUGE_VAL);
}

void ParticleVelsTest :: computeCacheForTest (void)
{
  
  Particle* p = new Particle();

  m_particleVels -> computeCacheFor(p);

  // looks a bit pointless to me that test...
  for(size_t i = 0; i < SPACE_DIMS; ++i)
    // that's what FakeFunctionParticle::computeCacheFor(p) should have
    // done:
    CPPUNIT_ASSERT_EQUAL ((double) i, p->v[i]);

}

