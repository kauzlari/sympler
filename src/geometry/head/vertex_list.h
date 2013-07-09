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



#ifndef __VERTEX_LIST_H
#define __VERTEX_LIST_H

#include <vector>

#include "geometric_primitives.h"


/*!
 * Stores a list of vertices
 */
class VertexList
{
 protected:
  /*!
   * The list of vertices
   */
  vector<point_t> m_vertices;

 public:
  /*!
   * Constructor
   */
  VertexList();

  /*!
   * Destructor
   */
  virtual ~VertexList();

  /*!
   * Add a new vertex
   * @param p Add vertex at this point
   */
  virtual int addVertex(const point_t &p) {
    m_vertices.push_back(p);
    return m_vertices.size()-1;
  }
	
  /*!
   * Add a new vertex
   * @param x x-coordinate of the vertex
   * @param y y-coordinate of the vertex
   * @param z z-coordinate of the vertex
   */
  virtual int addVertex(double x, double y, double z) {
    point_t v = { { {x, y, z} } };
    return addVertex(v);
  }

  /*!
   * Find vertex
   * @param vertex Find vertex at this point
   * @param eps Threshold for the difference of two vertices to be identified as one
   */
  virtual int findVertex(const point_t &vertex, double eps = 0);

  /*!
   * Create a new vertex only if it does not yet exist, otherwise return the existing one
   * @param vertex Find/create vertex at this point
   */
  virtual int needVertex(const point_t &vertex) {
    int v;
    v = findVertex(vertex, g_geom_eps);
    if (v == -1)
      v = addVertex(vertex);

    return v;
  }

  /*!
   * Create a new vertex only if it does not yet exist, otherwise return the existing one
   * @param x x-coordinate of the vertex
   * @param y y-coordinate of the vertex
   * @param z z-coordinate of the vertex
   */
  virtual int needVertex(double x, double y, double z) {
    point_t v = { { {x, y, z} } };
    return needVertex(v);
  }

  /*!
   * Return the list of vertices
   */
  virtual vector<point_t> &vertices() {
    return m_vertices;
  }
	
  /*!
   * Return vertex with index \a i
   * @param i Index of the vertex to return
   */
  virtual const point_t vertex(int i) const {
    return m_vertices[i];
  }

  /*!
   * Rotate all coordinates by \a by, i.e. in 3 dimensions by = 3 will change nothing
   * @param by Change the coordinates by this value
   */
  virtual void shiftCoords(int by);
};

#endif
