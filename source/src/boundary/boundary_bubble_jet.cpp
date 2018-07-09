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



#include "boundary_bubble_jet.h"

#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"

/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryBubbleJet> boundary_bubble_jet("BoundaryBubbleJet");


#define FRAME_OFF 0

#define FLOW_DIR 2
#define PERP_DIR1 0
#define PERP_DIR2 1


//---- Constructors/Destructor ----

BoundaryBubbleJet::BoundaryBubbleJet(Phase *phase): BoundaryArbitrary(phase)
{
	init();
}


BoundaryBubbleJet::~BoundaryBubbleJet()
{
}



//---- Methods ----

void BoundaryBubbleJet::init()
{
  m_properties.setClassName("BoundaryBubbleJet");

  m_properties.setDescription
    ("Simple diffusor type geometry. Cuboid inlet which expands to a cuboid reservoir.\n\n"
     "            inletHeight : inletHeight :    diffusorHeight    \n"
     "                                       --------------------- \n"
     "                                      |                     |\n" 
     "            --------------------------                      |\n"
     "           |            :             :                     |\n"
     "inletWidth |    inlet   :             :                     | diffusorWidth\n"
     "           |            :             :                     |\n"
     "            --------------------------                      |\n"
     "                                      |                     |\n"
     "                                       --------------------- \n"
     );

  INTPC
    (axisDir, m_axis_dir, -1,
     "Direction of the symmetry axis."
     );

  DOUBLEPC
    (reservoirWidth, m_reservoir_width, 0,
     "Width of the inlet (see diagram above).");
  DOUBLEPC
    (reservoirHeight, m_reservoir_height, 0,
     "Height of the inlet (see diagram above).");
  DOUBLEPC
    (inletWidth, m_inlet_width, 0,
     "Width of the inlet (see diagram above).");
  DOUBLEPC
    (inletHeight, m_inlet_height, 0,
     "Height of the inlet (see diagram above).");
  DOUBLEPC
    (diffusorWidth, m_diffusor_width, 0,
     "Width of the diffusor (see diagram above).");
  DOUBLEPC
    (diffusorHeight, m_diffusor_height, 0,
     "Height of the diffusor (see diagram above).");

  m_axis_dir = 0;
  m_reservoir_width = 10;
  m_reservoir_height = 5;
  m_inlet_height = 2;
  m_inlet_width = 1;
  m_inlet_height = 2;
  m_diffusor_width = 20;
  m_diffusor_height = 40;
}


	
#define ADDRECT(a, b, c, d)  {						\
    MSG_DEBUG("BoundaryBubbleJet::read", "Adding: " << a << ", " << b << ", " << c << ", " << d); \
    m_container->addWall(new WallTriangle(m_container, m_container->reflector(), a, b, c)); \
    m_container->addWall(new WallTriangle(m_container, m_container->reflector(), c, d, a)); } while(0)
	
	// Was a, b, c | b, c, d
	
#define ADDVERTEX(v, a, b, c) {			\
    point_t p;					\
    p.x = a+FRAME_OFF;				\
    p.y = b+FRAME_OFF;				\
    p.z = c+FRAME_OFF;				\
    v = m_container->addVertex(p); } while(0)


void BoundaryBubbleJet::setup()
{
  int v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;
  int v16, v17, v18, v19, v20, v21, v22, v23;
  double height1, height2;

  BoundaryArbitrary::setup();

  m_container = new WallContainer(m_reflectors[0]);

  height2 = m_reservoir_height;

  ADDVERTEX(v0, -m_reservoir_width/2, -m_reservoir_width/2, 0);
  ADDVERTEX(v1, -m_reservoir_width/2, +m_reservoir_width/2, 0);
  ADDVERTEX(v2, +m_reservoir_width/2, +m_reservoir_width/2, 0);
  ADDVERTEX(v3, +m_reservoir_width/2, -m_reservoir_width/2, 0);

  ADDVERTEX(v4, -m_reservoir_width/2, -m_reservoir_width/2, height2);
  ADDVERTEX(v5, -m_reservoir_width/2, +m_reservoir_width/2, height2);
  ADDVERTEX(v6, +m_reservoir_width/2, +m_reservoir_width/2, height2);
  ADDVERTEX(v7, +m_reservoir_width/2, -m_reservoir_width/2, height2);

  ADDVERTEX(v8, -m_inlet_width/2, -m_inlet_width/2, height2);
  ADDVERTEX(v9, -m_inlet_width/2, +m_inlet_width/2, height2);
  ADDVERTEX(v10, +m_inlet_width/2, +m_inlet_width/2, height2);
  ADDVERTEX(v11, +m_inlet_width/2, -m_inlet_width/2, height2);

  height1 = height2 + m_inlet_height;

  ADDVERTEX(v12, -m_inlet_width/2, -m_inlet_width/2, height1);
  ADDVERTEX(v13, -m_inlet_width/2, +m_inlet_width/2, height1);
  ADDVERTEX(v14, +m_inlet_width/2, +m_inlet_width/2, height1);
  ADDVERTEX(v15, +m_inlet_width/2, -m_inlet_width/2, height1);

  ADDVERTEX(v16, -m_diffusor_width/2, -m_diffusor_width/2, height1);
  ADDVERTEX(v17, -m_diffusor_width/2, +m_diffusor_width/2, height1);
  ADDVERTEX(v18, +m_diffusor_width/2, +m_diffusor_width/2, height1);
  ADDVERTEX(v19, +m_diffusor_width/2, -m_diffusor_width/2, height1);

  height2 = height1 + m_diffusor_height;

  ADDVERTEX(v20, -m_diffusor_width/2, -m_diffusor_width/2, height2);
  ADDVERTEX(v21, -m_diffusor_width/2, +m_diffusor_width/2, height2);
  ADDVERTEX(v22, +m_diffusor_width/2, +m_diffusor_width/2, height2);
  ADDVERTEX(v23, +m_diffusor_width/2, -m_diffusor_width/2, height2);


  /* Side walls */

  ADDRECT(v0, v4, v7, v3);
  ADDRECT(v3, v7, v6, v2);
  ADDRECT(v2, v6, v5, v1);
  ADDRECT(v1, v5, v4, v0);

  ADDRECT( v8, v12, v15, v11);
  ADDRECT(v11, v15, v14, v10);
  ADDRECT(v10, v14, v13,  v9);
  ADDRECT( v9, v13, v12,  v8);

  ADDRECT(v16, v20, v23, v19);
  ADDRECT(v19, v23, v22, v18);
  ADDRECT(v18, v22, v21, v17);
  ADDRECT(v17, v21, v20, v16);

  /* Other walls */

  ADDRECT( v8, v11,  v7,  v4);
  ADDRECT(v11, v10,  v6,  v7);
  ADDRECT(v10,  v9,  v5,  v6);
  ADDRECT( v9,  v8,  v4,  v5);


  ADDRECT(v12, v16, v19, v15);
  ADDRECT(v18, v14, v15, v19);
  ADDRECT(v18, v17, v13, v14);
  ADDRECT(v17, v16, v12, v13);

  //  ADDRECT(v15, v19, v16, v12);
  //  ADDRECT(v19, v15, v14, v18);
  //  ADDRECT(v14, v13, v17, v18);
  //  ADDRECT(v13, v12, v16, v17);

  ADDRECT(v0, v3, v2, v1);
  ADDRECT(v21, v22, v23, v20);

  if (m_axis_dir == 0)
    m_container->shiftCoords(1);
  else if (m_axis_dir == 1)
    m_container->shiftCoords(2);

  m_container->updateBoundingBox();
}



//---- Cell subdivision ----

void BoundaryBubbleJet::setup(Simulation* sim, ManagerCell *mgr)
{
  region_t *r;
//   point_t old_size;
  bool_point_t periodic = {{{ false, false, false }}};

  BoundaryArbitrary::setup(sim, mgr);

  MSG_DEBUG
    ("BoundaryCuboid::setup",
     "corner1 = " << m_container->boundingBox().corner1 << ", " <<
     "corner2 = " << m_container->boundingBox().corner2);

  MSG_DEBUG("BoundaryCuboid::setup", "Initializing cells.");
    
  r = mgr->cellSubdivide
    (sim->maxCutoff,
     m_container->boundingBox().corner1,
     m_container->boundingBox().corner2,
     periodic,
     1);
    
  mgr->cellSubdivisionFinished();

  mgr->assignContainer(m_container);

  m_container->toVTK(m_geometry_filename);
}

