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


#include "simulation.h"
#include "manager_cell.h"
#include "boundary_arbitrary.h"


//---- Constructors/Destructor ----

BoundaryArbitrary::BoundaryArbitrary(Phase *phase): Boundary(phase)
{
  init();
}


BoundaryArbitrary::~BoundaryArbitrary()
{
}



//---- Methods ----

void BoundaryArbitrary::setup()
{
	Boundary::setup();
}


bool BoundaryArbitrary::isInside(point_t point)
{
  return m_container->isInside(point);
}


bool BoundaryArbitrary::isInside(cuboid_t cuboid, const double& range)
{
  return m_container->isInside(cuboid, range);
}


bool BoundaryArbitrary::isInWallRange(const double& range, const point_t& point)
{
  return m_container->isInWallRange(range, point);
}


void BoundaryArbitrary::setup(Simulation *sim, ManagerCell *mgr)
{
  point_t old_size;
  
  m_proposedSize = m_container->boundingBox().size();
  old_size = m_proposedSize;
  
  MSG_DEBUG("BoundaryArbitrary::setup", "old_size = " << old_size);
  
  // next will loop over all PCs for adjustement of the box size
  Boundary::setup(sim, mgr);
  
  /* First we should care about stretching, then about the frame and moving of the hole box
     the special helpers of concrete boundaries may be strechted later. I hope this is no problem. */
  
  if (!(old_size == m_proposedSize)) {
    point_t factor;
    
    for (int i = 0; i < SPACE_DIMS; i++) {
      factor[i] = m_proposedSize[i] / old_size[i];
    }
    
    MSG_DEBUG("BoundaryArbitrary::setup", "old_size = (" << old_size.x << ", " <<
              old_size.y << ", " << old_size.z << ")");
    
    m_container->stretchBy(factor);
  }
  
  // adding a frame if wished...
  bool addOrigin = false;
  bool addEndPoint = false;
  
  cuboid_t new_bounds = m_container->boundingBox();
  
  for(int i = 0; i < SPACE_DIMS; i++)
    {
      if(m_frontFrame[i])
	{
	  addOrigin = true;
	  new_bounds.corner1[i] -= thickness();
	}
      if(m_endFrame[i])
	{
	  addEndPoint = true;
	  new_bounds.corner2[i] += thickness();
	}
    }
  
  if(addOrigin) 
    m_container->addVertex(new_bounds.corner1);  
  
  if(addEndPoint)
    m_container->addVertex(new_bounds.corner2);
  
  m_container->updateBoundingBox();
  
  // Yes, this gives you the same as m_container->boundingBox() would do
  MSG_DEBUG("BoundaryArbitrary::setup", "end: newBounds = " << new_bounds.corner1 << ", " << new_bounds.corner2);
  
}


void BoundaryArbitrary::init()
{
  m_properties.setClassName("BoundaryArbitrary");
  
  STRINGPC
    (geometryFileName, m_geometry_filename,
     "Geometry is converted and written to the VTK file given here.");
  
  m_geometry_filename = "geometry.vtk";    
}
