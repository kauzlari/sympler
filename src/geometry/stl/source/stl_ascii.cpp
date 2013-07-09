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

#include "stl_ascii.h"

#include "wall_triangle.h"


double STLAsciiParser::dist_eps = 1e-5;


/*--- STLAsciiParser ---*/

STLAsciiParser::STLAsciiParser(WallContainer *w): m_wall_container(w)
{
}


STLAsciiParser::~STLAsciiParser()
{
}


WallContainer *STLAsciiParser::read(istream *s)
{
  string key, str;
  bool immediate_eof = true;

  m_stream = s;

  readKey("solid");
  finishLine();

  key = readString();
  while (key != "endsolid" && !m_stream->eof()) {
    immediate_eof = false;

    if (key == "facet") {
      point_t normal;
      point_t vertices[3];
      int_point_t vindices;
      WallTriangle *w;

      readKey("normal");

      normal = readPoint();
      finishLine();
      readKey("outer");
      readKey("loop");
      finishLine();
 
      for (int i = 0; i < 3; i++) {
	readKey("vertex");
	vertices[i] = readPoint();
	finishLine();
      }

      readKey("endloop");
      readKey("endfacet");

      for (int i = 0; i < 3; i++) {
	int j;

	j = m_wall_container->findVertex(vertices[i], dist_eps);

	if (j < 0) {
	  j = m_wall_container->addVertex(vertices[i]);

   if(vertices[i].y > 2) MSG_DEBUG("STLAsciiParser::read", "New vertex with y>4 at: " << vertices[i]);
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

      //            MSG_DEBUG("STLAsciiParser::read", "nf = " << normal << ", nw = " << w->normal());

      //            MSG_DEBUG("STLAsciiParser::read", "WallTriangle: " << vindices);
    } else
      throw STLAsciiError("STLAsciiParser::read", "Keyword 'endsolid' or 'facet' expected.");

    key = readString();
  }

  if (immediate_eof)
    throw gError
      ("STLAsciiParser::read",
       "Reached the end of the file while reading the first key. "
       "Are you sure this is a STL? Perhaps it is a binary STL?");

  m_wall_container->updateBoundingBox();

  return m_wall_container;
}


WallContainer *STLAsciiParser::read(string filename)
{
  ifstream s;

  s.open(filename.c_str());

  if (!s.is_open())
    throw STLAsciiError("STLAsciiParser::read", "Error opening file '" + filename + "'.");

  read(&s);
  s.close();

  return m_wall_container;
}


string STLAsciiParser::readString()
{
  string s;

  (*m_stream) >> skipws >> s;

  return s;
}


void STLAsciiParser::readKey(string key)
{
  string s = readString();

  if (s != key)
  {
    if(key == "solid")
      throw STLAsciiError("STLAsciiParser::readKey", "Keyword 'solid' expected. Is the STL-file really in ASCII-format?");
    else
      throw STLAsciiError("STLAsciiParser::readKey", "Keyword '" + key + "' expected.");
  }
}

int STLAsciiParser::readInt()
{
  int i;

  (*m_stream) >> skipws >> i;

  return i;
}


double STLAsciiParser::readDouble()
{
  double d;

  (*m_stream) >> skipws >> d;

  return d;
}


point_t STLAsciiParser::readPoint()
{
  point_t p;

  p.x = readDouble();
  p.y = readDouble();
  p.z = readDouble();

  return p;
}


void STLAsciiParser::finishLine()
{
  while (m_stream->get() != '\n' && !m_stream->eof());
}
