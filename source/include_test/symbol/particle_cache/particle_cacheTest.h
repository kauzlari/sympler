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



#ifndef __PARTICLE_CACHE_TEST_H
#define __PARTICLE_CACHE_TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "fake_particle_cache.h"

using namespace std;

class ParticleCacheTest : public CPPUNIT_NS :: TestFixture
{

	  CPPUNIT_TEST_SUITE (ParticleCacheTest);

	  CPPUNIT_TEST (initTest);
	  CPPUNIT_TEST (setupTest);
	  CPPUNIT_TEST (checkOutputSymbolExistenceTest);
	  CPPUNIT_TEST (cleanSymbolTest);

	  CPPUNIT_TEST_SUITE_END ();


 public:
  
  /*!
   * Initialize objects
   */
  virtual void setUp (void);

  /*!
   * Helper allowing to instantiate different \a ParticleCache children for
   * different children of class \a ParticleCacheTest
   */
  virtual void setUpParticleCache (void);

  /*!
   * Delete objects
   */
  virtual void tearDown (void);

  
 protected:
  
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
   * Test \a ParticleCache::checkOutputSymbolExistence()
   */
  virtual void checkOutputSymbolExistenceTest (void);

  /*!
   * Test \a ParticleCache::cleanSymbol()
   */
  virtual void cleanSymbolTest (void);

  /*!
   * Helper for testing correctness of className in \a PropertyList of
   * \a ParticleCache
   */
  virtual void testPropClassName (void);

  /*!
   * Helper for testing correct initialisation of "species" attribute of
   * \a ParticleCache
   */
  virtual void testSpeciesAttr (void);

  /*!
   * Instance of class to be tested
   */
  ParticleCache *m_particleCache;

  /*!
   * Instance of \a Simulation used to setup \a m_particleCache
   */
  Simulation* m_simulation;
  
 private:

};

#endif
