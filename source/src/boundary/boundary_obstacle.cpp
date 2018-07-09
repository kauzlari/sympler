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



#include "boundary_obstacle.h"

#include "manager_cell.h"
#include "wall_triangle.h"

/* Register this Bounday with the factory. */
const Boundary_Register<BoundaryObstacle> boundary_obstacle("BoundaryObstacle");


#define FRAME_OFF 0

//---- Constructors/Destructor ----

BoundaryObstacle::BoundaryObstacle(Phase *phase): BoundaryCuboid(phase)
{
	init();
}


BoundaryObstacle::~BoundaryObstacle()
{
}



//---- Methods ----

void BoundaryObstacle::init()
{
  m_properties.setClassName("BoundaryObstacle");

  m_properties.setDescription
    ("Creates a simple channel and puts a cuboid obstacle in the middle.");

  for (int i = 0; i < 3; i++) {
    /* Register all properties */
    m_properties.addProperty
      ("obstacleCenter" + string(1, 'X'+i),
       PropertyList::DOUBLE, &m_obstacle_center[i],
       new PLCDoubleGreater(0),
       string(1, 'x'+i) + "-position of the obstacle.");
    m_properties.addProperty
      ("obstacleSize" + string(1, 'X'+i), 
       PropertyList::DOUBLE, &m_obstacle_size[i],
       new PLCDoubleGreater(0),
       string(1, 'x'+i) + "-component of the size of the obstacle.");
        
    /* Set default values */
    m_obstacle_center[i] = 1;
    m_obstacle_size[i] = 0.2;
  }
}


/*	
#define ADDRECT(a, b, c, d)  { \
	MSG_DEBUG("BoundaryObstacle::read", "Adding: " << a << ", " << b << ", " << c << ", " << d); \
    m_container->addWall(new WallTriangle(m_container, m_container->reflector(), a, b, c)); \
    m_container->addWall(new WallTriangle(m_container, m_container->reflector(), c, d, a)); } while(0) 
*/

#define ADDRECT(a, b, c, d, refl_name)  {				\
    /*MSG_DEBUG("BoundaryCuboid::read", "Adding: " << a << ", " << b << ", " << c << ", " << d);*/ \
    m_container->addWall(new WallTriangle(m_container, findReflector(refl_name), a, b, c)); \
    m_container->addWall(new WallTriangle(m_container, findReflector(refl_name), c, d, a)); } while(0)


    /*	m_container->addWall(new WallTriangle(m_container, m_container->reflector(), a, b, c)); \
    	m_container->addWall(new WallTriangle(m_container, m_container->reflector(), c, d, a)); } \
		 while(0) */
	
	// Was a, b, c | b, c, d
	
#define ADDVERTEX(v, a, b, c) {			\
    point_t p;					\
    p.x = a+FRAME_OFF;				\
    p.y = b+FRAME_OFF;				\
    p.z = c+FRAME_OFF;				\
    v = m_container->addVertex(p); } while(0)


void BoundaryObstacle::setup()
{
  int v0, v1, v2, v3, v4, v5, v6, v7;

  BoundaryCuboid::setup();

  // consistency check
  double upperObstcl, lowerObstcl;
  for(size_t i = 0; i < SPACE_DIMS; ++i)
    {
      upperObstcl = m_obstacle_center[i] + m_obstacle_size[i]/2;
      lowerObstcl = m_obstacle_center[i] - m_obstacle_size[i]/2;
      /*    MSG_DEBUG("BoundaryObstacle::setup", "upperObstcl=" << upperObstcl << ", lowerObstcl=" << lowerObstcl << ", m_box_size[i]=" << m_box_size[i]);*/
      if(upperObstcl > m_box_size[i] || lowerObstcl < 0)
	throw gError("BoundaryObstacle::setup", "Obstacle exceeds the box in " 
		     + string(1, 'x'+1) + "-direction");
    }
    
  ADDVERTEX(v0, 
	    m_obstacle_center.x - m_obstacle_size.x/2, 
	    m_obstacle_center.y - m_obstacle_size.y/2, 
	    m_obstacle_center.z - m_obstacle_size.z/2);
  ADDVERTEX(v1, 
	    m_obstacle_center.x + m_obstacle_size.x/2, 
	    m_obstacle_center.y - m_obstacle_size.y/2, 
	    m_obstacle_center.z - m_obstacle_size.z/2);
  ADDVERTEX(v2, 
	    m_obstacle_center.x + m_obstacle_size.x/2, 
	    m_obstacle_center.y + m_obstacle_size.y/2, 
	    m_obstacle_center.z - m_obstacle_size.z/2);
  ADDVERTEX(v3, 
	    m_obstacle_center.x - m_obstacle_size.x/2, 
	    m_obstacle_center.y + m_obstacle_size.y/2, 
	    m_obstacle_center.z - m_obstacle_size.z/2);

  ADDVERTEX(v4, 
	    m_obstacle_center.x - m_obstacle_size.x/2, 
	    m_obstacle_center.y - m_obstacle_size.y/2, 
	    m_obstacle_center.z + m_obstacle_size.z/2);
  ADDVERTEX(v5, 
	    m_obstacle_center.x + m_obstacle_size.x/2, 
	    m_obstacle_center.y - m_obstacle_size.y/2, 
	    m_obstacle_center.z + m_obstacle_size.z/2);
  ADDVERTEX(v6, 
	    m_obstacle_center.x + m_obstacle_size.x/2, 
	    m_obstacle_center.y + m_obstacle_size.y/2, 
	    m_obstacle_center.z + m_obstacle_size.z/2);
  ADDVERTEX(v7, 
	    m_obstacle_center.x - m_obstacle_size.x/2, 
	    m_obstacle_center.y + m_obstacle_size.y/2, 
	    m_obstacle_center.z + m_obstacle_size.z/2);

  ADDRECT(v2, v6, v5, v1, m_reflector_str[0]);
  ADDRECT(v4, v7, v3, v0, m_reflector_str[0]);

  ADDRECT(v1, v5, v4, v0, m_reflector_str[1]);
  ADDRECT(v7, v6, v2, v3, m_reflector_str[1]);

  ADDRECT(v3, v2, v1, v0, m_reflector_str[2]);
  ADDRECT(v5, v6, v7, v4, m_reflector_str[2]);

  m_container->updateBoundingBox();
	
  MSG_DEBUG("BoundaryObstacle::read", "box: " << m_container->boundingBox().corner1 << " - " << m_container->boundingBox().corner2);
}



