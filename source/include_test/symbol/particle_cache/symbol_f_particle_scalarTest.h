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


#ifndef __SYMBOL_F_PARTICLE_SCALAR_TEST_H
#define __SYMBOL_F_PARTICLE_SCALAR_TEST_H

#include "symbol_f_particle_scalar.h"

#include "symbol_f_particle_arbitraryTest.h"

using namespace std;


class SymbolFParticleScalarTest : public SymbolFParticleArbitraryTest
{
  CPPUNIT_TEST_SUITE (SymbolFParticleScalarTest);

  CPPUNIT_TEST (computeCacheForTest);
  CPPUNIT_TEST (setupTest);
  // next 5 defined in parent class only
  // FIXME: Does CPPUNIT allow for a better solution? I don't like it
  // that we do not see the definition of the tested functions in the
  // child class SymbolFParticleScalar (which is in principle a good
  // thing!), but we must(?) nonetheless
  // activate the tests here in the child class SymbolFParticleScalarTest
  // since the parent test class is an abstract class.
  CPPUNIT_TEST (initTest);
  CPPUNIT_TEST (setupOffsetTest);
  CPPUNIT_TEST (precomputeTest);
  CPPUNIT_TEST (registerWithParticleTest);
  CPPUNIT_TEST (setForceIndexToTest);

  CPPUNIT_TEST_SUITE_END ();
  
 public:

  /*!
   * Initialize objects
   */
  void setUp (void);

  
 protected:
  
  /*!
   * Test \a SymbolFParticleScalar::computeCacheFor(Particle* p)
   */
  void computeCacheForTest (void);

  /*!
   * Test \a SymbolFParticleScalar::init(Particle* p)
   */
  virtual void initTest (void); 

  /*!
   * Test \a SymbolFParticleScalar::setup()
   * FIXME: This test only checks what can go wrong on the 
   * \a SymbolFParticleScalar level and not responsibilities of parent classes.
   * This is a good thing! So fix this by extending the test-hierarchy 
   * to the parents and NOT by adding more test code here!
   */
  void setupTest (void);
  
  
 private:
  
};

#endif
