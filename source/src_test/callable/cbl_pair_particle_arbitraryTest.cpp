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

#include "cbl_pair_particle_arbitraryTest.h"

// commented out since this is an abstract class testing an abstract
// class
// CPPUNIT_TEST_SUITE_REGISTRATION (CblPairParticleArbitraryTest);

void CblPairParticleArbitraryTest :: initTest (void)
{
  /*!
   * Test init function. init() should have been called in constructor.
   */
  
  const PropertyList& properties = m_callable -> returnProperties();

  CPPUNIT_ASSERT_EQUAL (true, properties.exists("species1"));
  CPPUNIT_ASSERT_EQUAL (string("undefined"), m_callable -> returnSpeciesPair().first);
  CPPUNIT_ASSERT_EQUAL (string("undefined"), m_callable -> returnSpeciesPair().second);

  CPPUNIT_ASSERT_EQUAL (true, properties.exists("symbol"));
  CPPUNIT_ASSERT_EQUAL (string("undefined"), m_callable -> returnSymbolName());

  CPPUNIT_ASSERT_EQUAL (true, properties.exists("expression"));
  CPPUNIT_ASSERT_EQUAL (true, properties.exists("particleFactor_i"));
  CPPUNIT_ASSERT_EQUAL (true, properties.exists("particleFactor_j"));
  // no check of expression values since could be reset by init() of
  // child class  

  CPPUNIT_ASSERT_EQUAL (true, properties.exists("symmetry"));
  CPPUNIT_ASSERT_EQUAL ((size_t) 1, m_callable -> returnSymmetry());

  CPPUNIT_ASSERT_EQUAL (true, properties.exists("cutoff"));
  CPPUNIT_ASSERT_EQUAL (0., m_callable -> returnCutoff());

  CPPUNIT_ASSERT_EQUAL (true, properties.exists("overwrite"));
  CPPUNIT_ASSERT_EQUAL (false, m_callable -> returnOverwriteFlag());
  
}

