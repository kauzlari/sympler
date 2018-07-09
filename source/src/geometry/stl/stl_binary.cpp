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



#include <fstream>

#include "stl_binary.h"

#include "wall_triangle.h"


double STLBinaryParser::dist_eps = 1e-5;


/*--- STLBinaryParser ---*/

STLBinaryParser::STLBinaryParser(WallContainer *w): m_wall_container(w)
{
}


STLBinaryParser::~STLBinaryParser()
{
}


WallContainer *STLBinaryParser::read(istream *s)
{
  char str[100];
  size_t n_facets;
  math_vector_t<float> fnormal;
  math_vector_t<float> fvertices[3];
  unsigned short attr_size;

  /* Make sure the types we are using have the same size
     as the types in the STL specification. */

//   assert(sizeof(size_t) == 4); // 8 on 64 bit
  assert(sizeof(float) == 4);
  assert(sizeof(unsigned short) == 2);

  /* Read header */
  s->read(str, 80);
  str[80] = 0;
  MSG_DEBUG("STLBinaryParser::read", "STL header is " << endl << str);

  /* Number of facets */

  s->read((char*) &n_facets, 4/*not 64bit compatible: sizeof(size_t)*/);
  
  // on 64bit machines, the following removes rubbish in the upper unassigned 4Bytes
  // on 32bit machines, it's a bitwise "AND" with "1s" only 
  // I hope this is also BIG/LITTLE-Endian-independent
  size_t cleaner = 0xFFFFFFFF/*4294967295*/;
  /*the cleaner is 2^32-1; the decimal form produces a warning*/
  n_facets = n_facets&cleaner;
  
  MSG_DEBUG("STLBinaryParser::read", n_facets << " facets found.");

  for (size_t f = 0; f < n_facets; ++f) {
    int_point_t vindices;
    point_t normal, vertices[3];
    WallTriangle *w;

    s->read((char*) &fnormal.x, sizeof(float));
    s->read((char*) &fnormal.y, sizeof(float));
    s->read((char*) &fnormal.z, sizeof(float));

    s->read((char*) &fvertices[0].x, sizeof(float));
    s->read((char*) &fvertices[0].y, sizeof(float));
    s->read((char*) &fvertices[0].z, sizeof(float));

    s->read((char*) &fvertices[1].x, sizeof(float));
    s->read((char*) &fvertices[1].y, sizeof(float));
    s->read((char*) &fvertices[1].z, sizeof(float));

    s->read((char*) &fvertices[2].x, sizeof(float));
    s->read((char*) &fvertices[2].y, sizeof(float));
    s->read((char*) &fvertices[2].z, sizeof(float));

    s->read((char*) &attr_size, sizeof(unsigned short));

/*    // remove 64bit-rubbish again (see explanation above)
    fnormal.x = fnormal.x&4294967295;
    fnormal.y = fnormal.y&4294967295;
    fnormal.z = fnormal.z&4294967295;
    
    fvertices[0].x = (fvertices[0].x)&4294967295;
    fvertices[0].y = (fvertices[0].y)&4294967295;
    fvertices[0].z = (fvertices[0].z)&4294967295;
    
    fvertices[1].x = (fvertices[1].x)&4294967295;
    fvertices[1].y = (fvertices[1].y)&4294967295;
    fvertices[1].z = (fvertices[1].z)&4294967295;
    
    fvertices[2].x = (fvertices[2].x)&4294967295;
    fvertices[2].y = (fvertices[2].y)&4294967295;
    fvertices[2].z = (fvertices[2].z)&4294967295;*/
    
    /* Type conversion */
    for (size_t i = 0; i < SPACE_DIMS; ++i) {
      normal[i] = fnormal[i];
      for (size_t j = 0; j < 3; ++j)
        vertices[j][i] = fvertices[j][i];
    }

    for (size_t i = 0; i < 3; ++i) {
      int j;

      j = m_wall_container->findVertex(vertices[i], dist_eps);

      if (j < 0) {
        j = m_wall_container->addVertex(vertices[i]);
        //                    MSG_DEBUG("STLBinaryParser::read", "New vertex at: " << vertices[i]);
      }

      vindices[i] = j;
    }

    if (m_invert_normals)
      w = new WallTriangle
        (m_wall_container, m_wall_container->reflector(), vindices[2], vindices[1], vindices[0]);
    else
      w = new WallTriangle
        (m_wall_container, m_wall_container->reflector(), vindices[0], vindices[1], vindices[2]);

    m_wall_container->addWall(w);

    //            MSG_DEBUG("STLBinaryParser::read", "nf = " << normal << ", nw = " << w->normal());

    //            MSG_DEBUG("STLBinaryParser::read", "WallTriangle: " << vindices);
  } 

  m_wall_container->updateBoundingBox();
    
  return m_wall_container;
}


WallContainer *STLBinaryParser::read(string filename)
{
  ifstream s;

  s.open(filename.c_str());

  if (!s.is_open())
    throw STLBinaryError("STLBinaryParser::read", "Error opening file '" + filename + "'.");

  read(&s);
  s.close();

  return m_wall_container;
}


