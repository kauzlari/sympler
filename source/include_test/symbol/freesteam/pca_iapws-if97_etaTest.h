/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


#ifndef __PCA_IAPWSIF97_ETA_TEST_H
#define __PCA_IAPWSIF97_ETA_TEST_H

#ifdef HAVE_FREESTEAM

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "pca_iapws-if97_2varTest.h"
#include "pca_iapws-if97_eta.h"
#include "simulation.h"

using namespace std;


class PCacheIAPWSIF97etaTest : public PCacheIAPWSIF97TwoVarTest
{
  CPPUNIT_TEST_SUITE (PCacheIAPWSIF97etaTest);
  // implemented in parent class
  CPPUNIT_TEST (setupLUTTest);
  // implemented in parent class
  CPPUNIT_TEST (copyMySelfTest);
  CPPUNIT_TEST (calculateResultTest);
  CPPUNIT_TEST_SUITE_END ();  
  
 public:

  virtual void setUp (void);
  virtual void tearDown (void);

 protected:

  virtual void calculateResultTest (void);

 private:
  
};

#endif    // HAVE_FREESTEAM

#endif
