/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
 * David Kauzlaric <david.kauzlaric@frias.uni-freiburg.de>,
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



#include "smart_enumTest.h"
// #include <vector>

REGISTER_SMART_ENUM
(Int_Factory,
 "A factory for integers for testing the SmartEnum class for the int 'class'. See header file for factory class. Any test of SmartEnum for further types or classes must be done the same way."
);

CPPUNIT_TEST_SUITE_REGISTRATION (SmartEnumTest);

/*!
 * Initialize objects
 */
void SmartEnumTest :: setUp (void)
{
  smartEnumInt = new SmartEnum<Int_Factory>("TestForInt");
}

/*!
 * Delete objects
 */
void SmartEnumTest :: tearDown (void) 
{
  delete smartEnumInt;
}

void SmartEnumTest :: SmartEnumConstructorStringTest (void)
{
  /*!
   * Test constructor with string argument (constructor with empty argument list not defined)
   */
  CPPUNIT_ASSERT_EQUAL (string("TestForInt"), smartEnumInt->name());
  CPPUNIT_ASSERT_EQUAL (0, smartEnumInt->ordinal());
  CPPUNIT_ASSERT_EQUAL ((size_t) 1, smartEnumInt->mapF().size());

  // move this one also to setUp, if you need it also in another test
  SmartEnum<Int_Factory>* anotherSmartEnumInt = new SmartEnum<Int_Factory>("anotherTestForInt");

  CPPUNIT_ASSERT_EQUAL (string("anotherTestForInt"), anotherSmartEnumInt->name());
  CPPUNIT_ASSERT_EQUAL (1, anotherSmartEnumInt->ordinal());
  CPPUNIT_ASSERT_EQUAL (2, smartEnumInt->cardinality());
  CPPUNIT_ASSERT_EQUAL (2, anotherSmartEnumInt->cardinality());
  CPPUNIT_ASSERT_EQUAL ((size_t) 2, smartEnumInt->mapF().size());
  CPPUNIT_ASSERT_EQUAL ((size_t) 2, anotherSmartEnumInt->mapF().size());

  delete anotherSmartEnumInt;

}

