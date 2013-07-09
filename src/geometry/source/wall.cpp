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



#include "wall.h"
#include "wall_container.h"


//---- Constructors/Destructor ----

Wall::Wall(WallContainer *container, Reflector *reflector)
    : NodeOneChild((Node*) container), m_reflector(reflector)
{
}


Wall::~Wall()
{
}


Node *Wall::instantiateChild(const string &name)
{
	return NULL;
//    return Reflector_Factory::byName(name).instantiate(this);
}

void Wall::setup()
{
  NodeOneChild::setup();
/*  assert(m_parent);
  m_periodicity = ((WallContainer*) m_parent)->periodicityFront();*/
  const cuboid_t* tempBox = &(((WallContainer*) m_parent)->boundingBox());
  m_boxSize = tempBox->corner2 - tempBox->corner1;
//   MSG_DEBUG("Wall::setup", "corner2 = " << tempBox->corner2 << ", corner1 = " << tempBox->corner1 <<",  m_boxSize = " << m_boxSize);
}

void Wall::setBoundaryData(bool_point_t per, point_t size)
{
  m_boxSize = size;
  m_periodicity = per;
}
