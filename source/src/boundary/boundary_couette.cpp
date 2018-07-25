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



#include "boundary_couette.h"

#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"

/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryCouette> boundary_couette("BoundaryCouette");


//---- Constructors/Destructor ----

BoundaryCouette::BoundaryCouette(Phase *phase): BoundaryArbitrary(phase)
{
	init();
}


BoundaryCouette::~BoundaryCouette()
{
}

//---- Methods ----

void BoundaryCouette::init()
{
  m_properties.setClassName("BoundaryCouette");

  m_properties.setDescription("Creates a boundary for rotational Couette flow.");

  DOUBLEPC
    (innerRadius, m_inner_radius, 0,
     "Radius of the inner drum.");

  DOUBLEPC
    (outerRadius, m_outer_radius, 0,
     "Radius of the outer drum.");

  DOUBLEPC
    (height, m_height, 0,
     "Height of the two drum geometry.");

  BOOLPC
    (periodic, m_periodic,
     "Determines if the Couette drum should be continued periodically.");

  STRINGPC
    (innerReflector, m_inner_reflector,
     "Reflector to be used for particles hitting the inner drum.");

  STRINGPC
    (outerReflector, m_outer_reflector,
     "Reflector to be used for particles hitting the outer drum.");

  INTPC
    (nEdgeElements, m_n_edge_elements, 0,
     "Specifies how many (straight) elements to use in order to approximate "
     "the curvature of the inner and outer drum.");

  m_inner_radius = 1;
  m_outer_radius = 10;
  m_height = 10;
  m_periodic = true;
  m_inner_reflector = "default";
  m_outer_reflector = "default";
  m_n_edge_elements = 20;
}


	
#define ADDRECT(a, b, c, d, refl)  {                                    \
  /*MSG_DEBUG("BoundaryCouette::read", "Adding: " << a << ", " << b << ", " << c << ", " << d);*/ \
  m_container->addWall(new WallTriangle(m_container, refl, a, b, c));     \
  m_container->addWall(new WallTriangle(m_container, refl, c, d, a)); } while(0)


void BoundaryCouette::setup()
{
  bool_point_t periodic;
  Reflector *inner_refl, *outer_refl;

  MSG_DEBUG("BoundaryCouette::setup", "Creating geometry for rotational Couette flow.");

  BoundaryArbitrary::setup();

  if (!m_periodic)
    throw gError("BoundaryCouette::setup", "Non-periodic Couette flow is not yet supported.");

  periodic.x = periodic.y = false;
  periodic.z = m_periodic;

  /* Create a new container and tell it in which direction
     the boundary is periodic.
  */
	m_container = new WallContainer(m_reflectors[0]);
  m_container->setPeriodicityFront(periodic);
  m_container->setPeriodicityBack(periodic);

  inner_refl = findReflector(m_inner_reflector);
  outer_refl = findReflector(m_outer_reflector);

  /* Create the geometry. Loop over elements of the drums edges.
     In the first loop, the vertices are being created and in the
     second loop the vertices are connected by walls.
   */

  /* Create vertices. Note: The vertices are numbered in the order
     of their creation.
   */
  for (int i = 0; i < m_n_edge_elements; ++i) {
    point_t inner_pos, outer_pos;

    inner_pos.z = outer_pos.z = 0;

    inner_pos.x = m_inner_radius * cos(i*2*M_PI/m_n_edge_elements);
    inner_pos.y = m_inner_radius * sin(i*2*M_PI/m_n_edge_elements);
    outer_pos.x = m_outer_radius * cos(i*2*M_PI/m_n_edge_elements);
    outer_pos.y = m_outer_radius * sin(i*2*M_PI/m_n_edge_elements);

    m_container->addVertex(inner_pos);
    m_container->addVertex(outer_pos);

    inner_pos.z += m_height;
    outer_pos.z += m_height;

    m_container->addVertex(inner_pos);
    m_container->addVertex(outer_pos);
  }

  /* Create walls */
  for (int i = 0; i < m_n_edge_elements; ++i) {
    /* vib = Vertex Inner Bottom */
    int next_i, vib1, vob1, vit1, vot1, vib2, vob2, vit2, vot2;

    next_i = (i+1)%m_n_edge_elements;

    vib1 = i*4;
    vob1 = i*4+1;
    vit1 = i*4+2;
    vot1 = i*4+3;

    vib2 = next_i*4;
    vob2 = next_i*4+1;
    vit2 = next_i*4+2;
    vot2 = next_i*4+3;

    ADDRECT(vib1, vib2, vit2, vit1, inner_refl);
    ADDRECT(vot1, vot2, vob2, vob1, outer_refl);
  }

	m_container->updateBoundingBox();
}



//---- Cell subdivision ----

void BoundaryCouette::setup(Simulation* sim, ManagerCell *mgr)
{
  region_t *r;
  bool_point_t periodic;

  periodic.x = periodic.y = false;
  periodic.z = m_periodic;

  BoundaryArbitrary::setup(sim, mgr);
    
  MSG_DEBUG
    ("BoundaryCouette::setup",
     "corner1 = " << m_container->boundingBox().corner1 << ", " <<
     "corner2 = " << m_container->boundingBox().corner2);

  MSG_DEBUG("BoundaryCouette::setup", "Initializing cells.");
    
  r = mgr->cellSubdivide
    (sim->maxCutoff, m_container->boundingBox().corner1, m_container->boundingBox().corner2, periodic);
    
  mgr->cellSubdivisionFinished();

  mgr->assignContainer(m_container);

  m_container->toVTK(m_geometry_filename);
}
