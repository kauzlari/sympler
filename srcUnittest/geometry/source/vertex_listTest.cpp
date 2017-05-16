/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2013, 
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



#include "vertex_listTest.h"
#include <vector>

CPPUNIT_TEST_SUITE_REGISTRATION (vertex_listTest);

void vertex_listTest :: setUp (void)
{
  /*!
   * Initialize objects
   */
  a = new VertexList();
  b = new VertexList();
}

void vertex_listTest :: tearDown (void) 
{
  /*!
   * Delete objects
   */
  delete a; delete b;
}

void vertex_listTest :: addVertexTest (void)
{
  /*!
   * Test adding vertices and returning its index
   */
  point_t v = { { {1, 2, 3} } };
  point_t w = { { {4, 5, 6} } };
  CPPUNIT_ASSERT_EQUAL (0, a->addVertex(v));
  CPPUNIT_ASSERT_EQUAL (1, a->addVertex(w));
  CPPUNIT_ASSERT_EQUAL (2, a->addVertex(0, 0, 0));
}

void vertex_listTest :: findVertexTest (void)
{
  /*!
   * Test finding vertices with eps = 0 and eps != 0 (distance)
   */
  point_t v = { { {1, 2, 3} } };
  point_t w = { { {4, 5, 6} } };
  point_t x = { { {2, 2, 3} } };
  a->addVertex(v);
  a->addVertex(w);
  CPPUNIT_ASSERT_EQUAL (0, a->findVertex(v, 0));
  CPPUNIT_ASSERT_EQUAL (-1, b->findVertex(v, 0));
  CPPUNIT_ASSERT_EQUAL (1, a->findVertex(w, 0));
  CPPUNIT_ASSERT_EQUAL (-1, a->findVertex(x, 0));

  CPPUNIT_ASSERT_EQUAL (-1, a->findVertex(x, 0.5));
  CPPUNIT_ASSERT_EQUAL (0, a->findVertex(x, 1.5));
}

void vertex_listTest :: needVertexTest (void)
{
  /*!
   * Test creating a new vertex only if it does not yet exist, otherwise returning the existing one
   */
  point_t v = { { {1, 2, 3} } };
  point_t w = { { {4, 5, 6} } };
  a->addVertex(v);
  a->addVertex(w);
  CPPUNIT_ASSERT_EQUAL (0, a->needVertex(v));
  CPPUNIT_ASSERT_EQUAL (1, a->needVertex(w));
  CPPUNIT_ASSERT_EQUAL (2, a->needVertex(7, 8, 9));
}

void vertex_listTest :: shiftCoordsTest (void)
{
  /*!
   * Test rotating coordinates
   */
  point_t v = { { {1, 2, 3} } };
  point_t w = { { {2, 3, 1} } };
  point_t x = { { {3, 1, 2} } };
  point_t vv = { { {4, 5, 6} } };
  point_t ww = { { {5, 6, 4} } };
  point_t xx = { { {6, 4, 5} } };
  a->addVertex(v);
  a->addVertex(vv);
  a->shiftCoords(0);
  CPPUNIT_ASSERT_EQUAL (v, a->vertex(0));
  CPPUNIT_ASSERT_EQUAL (vv, a->vertex(1));
  a->shiftCoords(1);
  CPPUNIT_ASSERT_EQUAL (w, a->vertex(0));
  CPPUNIT_ASSERT_EQUAL (ww, a->vertex(1));
  a->shiftCoords(1);
  CPPUNIT_ASSERT_EQUAL (x, a->vertex(0));
  CPPUNIT_ASSERT_EQUAL (xx, a->vertex(1));
  a->shiftCoords(1);
  CPPUNIT_ASSERT_EQUAL (v, a->vertex(0));
  CPPUNIT_ASSERT_EQUAL (vv, a->vertex(1));
  a->shiftCoords(3);
  CPPUNIT_ASSERT_EQUAL (v, a->vertex(0));
  CPPUNIT_ASSERT_EQUAL (vv, a->vertex(1));
  a->shiftCoords(9);
  CPPUNIT_ASSERT_EQUAL (v, a->vertex(0));
  CPPUNIT_ASSERT_EQUAL (vv, a->vertex(1));
}


