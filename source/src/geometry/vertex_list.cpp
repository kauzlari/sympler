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



#include <cassert>

#include "misc.h"
#include "vertex_list.h"

/*--- VertexList ---*/

VertexList::VertexList()
{
}

VertexList::~VertexList()
{
}

int VertexList::findVertex(const point_t &vertex, double eps)
{
    int index = -1;

    if (eps == 0) {
        for (size_t i = 0; i < m_vertices.size(); i++) {
            if (m_vertices[i] == vertex)
            {
              index = i;
              assert(index > -1);
            }
        }
    } else {
        for (size_t i = 0; i < m_vertices.size(); i++) {
            if ((m_vertices[i] - vertex).abs() < eps)
            {
              index = i;
              assert(index > -1);
            }
        }
    }

    return index;
}


void VertexList::shiftCoords(int by)
{
  for (vector<point_t>::iterator i = m_vertices.begin(); i != m_vertices.end(); ++i) {
    point_t help = *i;

    for (int j = 0; j < SPACE_DIMS; ++j)
      (*i)[j] = help[(j+by)%SPACE_DIMS];
  }
}
