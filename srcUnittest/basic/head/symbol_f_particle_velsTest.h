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



#ifndef __SYMBOL_F_PARTICLE_VELS_TEST_H
#define __SYMBOL_F_PARTICLE_VELS_TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "symbol_f_particle_vels.h"

using namespace std;

class SymbolFParticleVelsTest : public CPPUNIT_NS :: TestFixture
{
  CPPUNIT_TEST_SUITE (SymbolFParticleVelsTest);
  CPPUNIT_TEST (computeCacheForTest);
  CPPUNIT_TEST (initTest);
  CPPUNIT_TEST (setupTest);
  CPPUNIT_TEST (setupOffsetTest);
  CPPUNIT_TEST_SUITE_END ();
  
 public:
  /*!
   * Initialize objects
   */
  void setUp (void);

  /*!
   * Delete objects
   */
  void tearDown (void);

  
 protected:
  
  /*!
   * Test \a SymbolFParticleVels::computeCacheFor(Particle* p)
   */
  void computeCacheForTest (void);
  
  /*!
   * Test \a SymbolFParticleVels::init(Particle* p)
   */
  void initTest (void);
  
  /*!
   * Test \a SymbolFParticleVels::setup()
   * FIXME: This test only checks what can go wrong on the 
   * \a SymbolFParticleVels level and not responsibilities of parent classes.
   * This is a good thing! So fix this by extending the test-hierarchy 
   * to the parents and NOT by adding more test code here!
   */
  void setupTest (void);
  
  /*!
   * Test \a SymbolFParticleVels::setupOffset(Particle* p)
   */
  void setupOffsetTest (void);
  
 private:
  
  SymbolFParticleVels *m_symbolFParticleVels;

};

#endif
