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



#ifndef __SYMBOL_F_PARTICLE_ARBITRARY_TEST_H
#define __SYMBOL_F_PARTICLE_ARBITRARY_TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "symbol_f_particle_arbitrary.h"

using namespace std;

class SymbolFParticleArbitraryTest : public CPPUNIT_NS :: TestFixture
{
  
 public:
  
  /*!
   * Initialize objects
   */
  virtual void setUp (void);

  /*!
   * Delete objects
   */
  virtual void tearDown (void);

  
 protected:

  /*!
   * Helper function; not a test; meant to be called by setUp() of 
   * child test class
   */
  void setUpFunctionParticle();
  
  /*!
   * Test \a SymbolFParticleArbitrary::init(Particle* p)
   */
  virtual void initTest (void);
  
  /*!
   * Test \a SymbolFParticleArbitrary::setup()
   * FIXME: This test only checks what can go wrong on the 
   * \a SymbolFParticleArbitrary level and not responsibilities of parent classes.
   * This is a good thing! So fix this by extending the test-hierarchy 
   * to the parents and NOT by adding more test code here!
   */
  virtual void setupTest (void);

  /*!
   * Test \a SymbolFParticleArbitrary::setupOffset(Particle* p)
   */
  virtual void setupOffsetTest (void);
  
  /*!
   * Test \a SymbolFParticleArbitrary::precompute()
   */
  virtual void precomputeTest (void);
    
  /*!
   * Test \a SymbolFParticleArbitrary::computeCacheFor(Particle* p)
   */
  virtual void computeCacheForTest (void) = 0;

  /*!
   * Helper for \a computeCacheForTest (void)
   */
  virtual void computeCacheFor(Particle* p, size_t forceIndex);
  
  /*!
   * Test \a SymbolFParticleArbitrary::registerWithParticle()
   * Since the function does not do anything, we do not test anything
   * so far. If a child of SymbolFParticleArbitrary calls the function
   * of the parent, then the test function of the test child should 
   * call this function here of its test parent
   */
  virtual void registerWithParticleTest (void) {
    m_symbolFParticleArbitrary -> registerWithParticle();
    // no test required so far since function does nothing
  }

  /*!
   * Test \a SymbolFParticleArbitrary::setForceIndexTo(size_t forceIndex)
   */
  virtual void setForceIndexToTest (void);

  /*!
   * Instance of class to be tested
   */
  SymbolFParticleArbitrary *m_symbolFParticleArbitrary;

  /*!
   * Instance of \a Simulation used to setup \a m_symbolFParticleArbitrary
   */
  Simulation* m_simulation;
  
 private:

};

#endif
