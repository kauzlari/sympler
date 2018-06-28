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



#ifndef __CBL_PAIR_PARTICLE_ARBITRARY_TEST_H
#define __CBL_PAIR_PARTICLE_ARBITRARY_TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "cbl_pair_particle_arbitrary.h"

using namespace std;


/*!
 * Abstract parent for test classes of children of abstract class 
 * \a CblPairParticleArbitrary
 * FIXME: setup() not tested yet!
 */
class CblPairParticleArbitraryTest : public CPPUNIT_NS :: TestFixture
{
  /* CPPUNIT_TEST_SUITE (CblPairParticleArbitraryTest); */
  /* CPPUNIT_TEST (initTest); */
  /* CPPUNIT_TEST_SUITE_END (); */
  
 public:
  
  /*!
   * Initialise used objects
   */
  virtual void setUp (void) = 0;
  
  /*!
   * Delete used objects
   */
  virtual void tearDown (void) = 0;
  
 protected:
  
  /*!
   * Test init function
   */
  virtual void initTest (void);

  /*!
   * Instance of tested class
   */
  CblPairParticleArbitrary *m_callable;
  
 private:
  
};

#endif
