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



#include "boundary_cuboid.h"

#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"

/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryCuboid> boundary_cuboid("BoundaryCuboid");


#define FRAME_OFF 0

//---- Constructors/Destructor ----

BoundaryCuboid::BoundaryCuboid(Phase *phase): BoundaryArbitrary(phase)
{
	init();
}


BoundaryCuboid::~BoundaryCuboid()
{
}

//---- Methods ----

void BoundaryCuboid::init()
{
  m_properties.setClassName("BoundaryCuboid");

  m_properties.setDescription("Defines a simple cuboid simulation region.");

  for (int i = 0; i < 3; i++) {
    /* Register all properties */
    m_properties.addProperty
      ("box" + string(1, 'X'+i), PropertyList::DOUBLE, 
       &m_box_size[i],
       new PLCDoubleGreater(0),
       string(1, 'x'+i) + "-size of the simulation box.");
    m_properties.addProperty
      ("periodic" + string(1, 'X'+i), PropertyList::BOOLEAN, &m_periodic[i], NULL,
       "Simulation box is periodic in " + string(1, 'x'+i) + "-direction.");

    m_properties.addProperty
      ("reflector" + string(1, 'X'+i), PropertyList::STRING, &m_reflector_str[i], NULL,
       "Name of the reflector for walls in " + string(1, 'x'+i) + "-direction.");
        
    /* Set default values */
    m_box_size[i] = 10;
    m_periodic[i] = true;
    m_reflector_str[i] = "default";
  }
}


	
#define ADDRECT(a, b, c, d, refl_name)  {                                         \
	/*MSG_DEBUG("BoundaryCuboid::read", "Adding: " << a << ", " << b << ", " << c << ", " << d);*/ \
	m_container->addWall(new WallTriangle(m_container, findReflector(refl_name), a, b, c)); \
	m_container->addWall(new WallTriangle(m_container, findReflector(refl_name), c, d, a)); } while(0)
	
	// Was a, b, c | b, c, d
	
#define ADDVERTEX(v, a, b, c) { \
	point_t p; \
	p.x = a+FRAME_OFF; \
	p.y = b+FRAME_OFF; \
	p.z = c+FRAME_OFF; \
	v = m_container->addVertex(p); } while(0)


void BoundaryCuboid::setup()
{
  int v0, v1, v2, v3, v4, v5, v6, v7;

//   MSG_DEBUG("BoundaryCuboid::setup", "...");

  BoundaryArbitrary::setup();

	m_container = new WallContainer(m_reflectors[0]);
  m_container->setPeriodicityFront(m_periodic);
  m_container->setPeriodicityBack(m_periodic);
		
  ADDVERTEX(v0, 0, 0, 0);
  ADDVERTEX(v1, m_box_size.x, 0, 0);
  ADDVERTEX(v2, m_box_size.x, m_box_size.y, 0);
  ADDVERTEX(v3, 0, m_box_size.y, 0);

  ADDVERTEX(v4, 0, 0, m_box_size.z);
  ADDVERTEX(v5, m_box_size.x, 0, m_box_size.z);
  ADDVERTEX(v6, m_box_size.x, m_box_size.y, m_box_size.z);
  ADDVERTEX(v7, 0, m_box_size.y, m_box_size.z);

  if (!m_periodic.x) {
    ADDRECT(v1, v5, v6, v2, m_reflector_str[0]);
    ADDRECT(v0, v3, v7, v4, m_reflector_str[0]);
  }
  if (!m_periodic.y) {
    ADDRECT(v0, v4, v5, v1, m_reflector_str[1]);
    ADDRECT(v3, v2, v6, v7, m_reflector_str[1]);
  }
  if (!m_periodic.z) {
    ADDRECT(v0, v1, v2, v3, m_reflector_str[2]);
    ADDRECT(v4, v7, v6, v5, m_reflector_str[2]);
  }

	m_container->updateBoundingBox();
}



//---- Cell subdivision ----

#define OUTLET 0
#define NEIGHBOR 1

void BoundaryCuboid::setup(Simulation* sim, ManagerCell *mgr)
{
  region_t *r;

  BoundaryArbitrary::setup(sim, mgr);
    
  MSG_DEBUG
    ("BoundaryCuboid::setup",
     "corner1 = " << m_container->boundingBox().corner1 << ", " <<
     "corner2 = " << m_container->boundingBox().corner2);

  MSG_DEBUG("BoundaryCuboid::setup", "Initializing cells.");
    
  r = mgr->cellSubdivide
    (sim->maxCutoff, m_container->boundingBox().corner1, m_container->boundingBox().corner2, m_periodic, 1);
    
  mgr->cellSubdivisionFinished();

  mgr->assignContainer(m_container);

  m_container->toVTK(m_geometry_filename);
}
