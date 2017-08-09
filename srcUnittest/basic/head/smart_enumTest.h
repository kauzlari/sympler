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



#ifndef SMART_ENUMTEST_H
#define SMART_ENUMTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "smart_enum.h"

using namespace std;


//---- Factories --- Done the very same way as in the actual code, for each enumed "class" we wish to test ----

/*!
 * Abstract factory class for integers
 */
class Int_Factory : public SmartEnum<Int_Factory> {
public:
	virtual int *instantiate() const = 0;

protected:
	Int_Factory(const string &name) :
		SmartEnum<Int_Factory>(name) {
	}
};

/*!
 * Actual registration class for classes derived from int. (Only intended for testing!)
 */
template <class T> class Int_Register : public Int_Factory {
public:
	Int_Register(const string &name) :
		Int_Factory(name) {
	}

	virtual int *instantiate() const;
};


class SmartEnumTest : public CPPUNIT_NS :: TestFixture
{
  CPPUNIT_TEST_SUITE (SmartEnumTest);
  CPPUNIT_TEST (SmartEnumConstructorStringTest);
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
     * Test of constructor with string argument
     */
    void SmartEnumConstructorStringTest (void);

  private:
    SmartEnum<Int_Factory> *smartEnumInt;
//    SmartEnum *b;
};


#endif
