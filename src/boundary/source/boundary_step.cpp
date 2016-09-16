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



#include "boundary_step.h"

#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"

/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryStep> boundary_step("BoundaryStep");


#define FRAME_OFF 0


//---- Constructors/Destructor ----

BoundaryStep::BoundaryStep(Phase *phase): BoundaryArbitrary(phase)
{
	init();
}


BoundaryStep::~BoundaryStep()
{
}



//---- Methods ----

void BoundaryStep::init()
{
	m_properties.setClassName("BoundaryStep");

  m_properties.setDescription(
    "Creates a simple step boundary.\n\n"
    "          : bigLflow : beforeStep : afterStep : smallLflow    \n"
    "           ------------------------------------------------             \n"
    "          |          :            :           :            |            \n"
    "          |          :            :           :   outlet   | smallHeight\n"
    "bigHeight |  inlet   :            :           :            |            \n"
    "          |          :           --------------------------             \n"
    "          |          :          |\n"
    "          |          :          |\n"
    "           --------------------- \n"
  );

  m_properties.addProperty
    ("sideWall", PropertyList::STRING, &m_side_wall,
     new PLCStringList2("periodic", "wall"),
     "Determines whether the direction perpendicular to the flow and perpedicular to the step "
     "has a wall on both sides or is treated periodically.");
  m_properties.addProperty
    ("flowDir", PropertyList::INT, &m_flow_dir,
     new PLCIntList3(0, 1, 2),
     "Direction of flow (0,_1_or_2).");
  m_properties.addProperty
    ("stepLoc", PropertyList::INT, &m_step_loc,
     new PLCIntList3(0, 1, 2),
     "Direction of the step (0,_1_or_2).");
  DOUBLEPC
    (bigHeight, m_big_height, 0,
     "bigHeight in the above diagram");
  DOUBLEPC
    (bigLflow, m_big_lflow, 0,
     "bigLflow in the above diagram");
  DOUBLEPC
    (beforeStep, m_before_step, 0,
     "beforeStep in the above diagram");
  DOUBLEPC
    (afterStep, m_after_step, 0,
     "afterStep in the above diagram");
  DOUBLEPC
    (smallLflow, m_small_lflow, 0,
     "smallLflow in the above diagram");
  DOUBLEPC
    (width, m_width, 0,
     "Width of the simulation region.");
  DOUBLEPC
    (smallHeight, m_small_height, 0,
     "smallHeight in the above diagram");
  DOUBLEPC
    (temperatureSides, m_temperature_sides, 0,
     "Obsolete.");
  DOUBLEPC
    (temperatureStepDir, m_temperature_step_dir, 0,
     "Obsolete.");

  m_side_wall = "periodic";
  m_flow_dir = 0;
  m_step_loc = 1;
  m_big_height = m_big_lflow = m_before_step = m_after_step = m_small_lflow = m_width = m_small_height = 0;
  m_temperature_sides = m_temperature_step_dir = 1;
}



#define ADDRECT(a, b, c, d)  { \
	/*MSG_DEBUG("BoundaryStep::read", "Adding: " << a << ", " << b << ", " << c << ", " << d);*/ \
	m_container->addWall(new WallTriangle(m_container, m_container->reflector(), a, b, c)); \
	m_container->addWall(new WallTriangle(m_container, m_container->reflector(), c, d, a)); } while(0)
	
	// Was a, b, c | b, c, d
	
#define ADDVERTEX(v, a, b, c) { \
	point_t p; \
	p.x = a+FRAME_OFF; \
	p.y = b+FRAME_OFF; \
	p.z = c+FRAME_OFF; \
	v = m_container->addVertex(p); } while(0)


void BoundaryStep::setup()
{
	int v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13;
  bool_point_t periodic;

	BoundaryArbitrary::setup();

  m_side_periodic = (m_side_wall == "periodic");

//	point_t o, s1, s2;
	
	m_third_dir = 0;
	
	while (m_third_dir == m_flow_dir || m_third_dir == m_step_loc)
		m_third_dir++;
		
	MSG_DEBUG("BoundaryStep::read", "m_third_dir = " << m_third_dir);

	m_container = new WallContainer(m_reflectors[0]);

  periodic[m_flow_dir] = true;
  periodic[m_step_loc] = false;
  periodic[m_third_dir] = m_side_periodic;
  m_container->setPeriodicityFront(periodic);
  m_container->setPeriodicityBack(periodic);

	ADDVERTEX(v0, 0, 0, 0);
	ADDVERTEX(v1, 0, m_big_height, 0);
	ADDVERTEX(v2, m_before_step+m_big_lflow, m_big_height, 0);
	ADDVERTEX(v3, m_before_step+m_small_lflow+m_big_lflow+m_after_step, m_big_height, 0);
	ADDVERTEX(v4, m_before_step+m_small_lflow+m_big_lflow+m_after_step, m_big_height-m_small_height, 0);
	ADDVERTEX(v5, m_before_step+m_big_lflow, m_big_height-m_small_height, 0);
	ADDVERTEX(v6, m_before_step+m_big_lflow, 0, 0);

	ADDVERTEX(v7, 0, 0, m_width);
	ADDVERTEX(v8, 0, m_big_height, m_width);
	ADDVERTEX(v9, m_before_step+m_big_lflow, m_big_height, m_width);
	ADDVERTEX(v10, m_before_step+m_small_lflow+m_big_lflow+m_after_step, m_big_height, m_width);
	ADDVERTEX(v11, m_before_step+m_small_lflow+m_big_lflow+m_after_step, m_big_height-m_small_height, m_width);
	ADDVERTEX(v12, m_before_step+m_big_lflow, m_big_height-m_small_height, m_width);
	ADDVERTEX(v13, m_before_step+m_big_lflow, 0, m_width);
	
  if (!m_side_periodic) {
    ADDRECT(v0, v6, v2, v1);
    ADDRECT(v5, v4, v3, v2);
    ADDRECT(v8, v9, v13, v7);
    ADDRECT(v9, v10, v11, v12);
  }
		
	ADDRECT(v7, v0, v1, v8);
	ADDRECT(v6, v13, v12, v5);
	ADDRECT(v4, v11, v10, v3);
	
	ADDRECT(v0, v7, v13, v6);
	ADDRECT(v5, v12, v11, v4);
	ADDRECT(v1, v3, v10, v8);

	m_container->updateBoundingBox();
	
	MSG_DEBUG("BoudaryStepStoch::read", "box: " << m_container->boundingBox().corner1 << " - " << m_container->boundingBox().corner2);
}



//---- Cell subdivision ----

void BoundaryStep::setup(Simulation* sim, ManagerCell *mgr)
{
  point_t old_size;

	old_size = m_container->boundingBox().corner2;
  /*	old_size = m_size;

	MSG_DEBUG("BoundaryStep::setup", "old_size = " << old_size);

	m_particle_creator->adjustBoxSize(m_size);
	m_volume = m_size.x * m_size.y * m_size.z;

	if (!(old_size == m_size)) {
  point_t factor;
	
  for (int i = 0; i < SPACE_DIMS; i++) {
  factor[i] = m_size[i] / old_size[i];
  }
		
  m_container->stretchBy(factor);
	
  m_big_height *= factor[m_step_loc];
  m_big_lflow *= factor[m_flow_dir];
  m_before_step *= factor[m_flow_dir]; 
  m_after_step *= factor[m_flow_dir];
  m_small_lflow *= factor[m_flow_dir];
  m_width *= factor[m_third_dir];
  m_small_height *= factor[m_step_loc];
  }*/

  BoundaryArbitrary::setup(sim, mgr);

	if (!(old_size == m_proposedSize)) {
		point_t factor;
	
		for (int i = 0; i < SPACE_DIMS; i++) {
			factor[i] = m_proposedSize[i] / old_size[i];
		}
		
		m_big_height *= factor[m_step_loc];
		m_big_lflow *= factor[m_flow_dir];
		m_before_step *= factor[m_flow_dir]; 
		m_after_step *= factor[m_flow_dir];
		m_small_lflow *= factor[m_flow_dir];
		m_width *= factor[m_third_dir];
		m_small_height *= factor[m_step_loc];
  }
	

/* Cell subdivision */
	
  point_t c1, c2;
  region_t *inlet, *middle, *outlet;
  bool_point_t p;

  p.x = p.y = p.z = false;
  p[m_third_dir] = m_side_periodic;
	
	MSG_DEBUG("BoundaryStep::setup", "Initalizing cells.");
    
  for (int i = 0; i < SPACE_DIMS; i++) {
    c1[i] = 0;
  }
    		
  Phase* phase = (Phase*) m_parent;
						
  c2[m_flow_dir] = m_big_lflow;
  // there are definitely no frozen particles in flow direction
  if(m_frontFrame[m_flow_dir]) throw gError("BoundaryStep::setup: cannot have wall particles" 
                                            " in flow direction");
  if(m_endFrame[m_flow_dir]) throw gError("BoundaryStep::setup: cannot have wall particles" 
                                          " in flow direction");
  c2[m_step_loc] = m_big_height;
  // next two lines necessary if frozen particles frame
  if(m_frontFrame[m_step_loc]) c1[m_step_loc] -= phase -> pairCreator() -> interactionCutoff();
  if(m_endFrame[m_step_loc]) c2[m_step_loc] += phase -> pairCreator() -> interactionCutoff();
//   if(m_frontFrame[m_step_loc]) c1[m_step_loc] -= sim -> maxCutoff;
//   if(m_endFrame[m_step_loc]) c2[m_step_loc] += sim -> maxCutoff;
  c2[m_third_dir] = m_width;
  // next two lines necessary if frozen particles frame
  if(m_frontFrame[m_third_dir])
    if(m_side_periodic) throw gError("BoundaryStep::setup: frontFrame for periodic side"
                                     " not possible");
    else c1[m_third_dir] -= phase -> pairCreator() -> interactionCutoff();
//     else c1[m_third_dir] -= sim -> maxCutoff;
  if(m_endFrame[m_third_dir])
    if(m_side_periodic) throw gError("BoundaryStep::setup: endFrame for periodic side"
                                     " not possible");
    else c2[m_third_dir] += phase -> pairCreator() -> interactionCutoff();
//     else c2[m_third_dir] += sim -> maxCutoff;

    
  p[m_flow_dir] = true;
  inlet = mgr->cellSubdivide(sim->maxCutoff, c1, c2, p, 1);
    
  c1[m_flow_dir] = c2[m_flow_dir];
  c2[m_flow_dir] = c2[m_flow_dir] + m_before_step + m_after_step;
    
  p[m_flow_dir] = false;
  middle = mgr->cellSubdivide(sim->maxCutoff, c1, c2, p, 0);
    
  c1[m_flow_dir] = c2[m_flow_dir];
  c2[m_flow_dir] = c2[m_flow_dir] + m_small_lflow;

  p[m_flow_dir] = true;
  outlet = mgr->cellSubdivide(sim->maxCutoff, c1, c2, p, 2);

  MSG_DEBUG("BoundaryStep::setup", "inlet->n_cells = " << inlet->n_cells << ", middle->n_cells = " << middle->n_cells << ", outlet->n_cells = " << outlet->n_cells);


  MSG_DEBUG("BoundaryStep::setup", "Establishing connections between inlet, middle part and outlet.");

  /* Remove periodicity from middle */
  //    ManagerCell::removePeriodicity(m_flow_dir, middle);
            
  /* Connect inlet to middle */
  mgr->connect(m_flow_dir, inlet, middle, p, OUTLET);
            
  /* Connect outlet to middle */
  mgr->connect(-m_flow_dir-1, outlet, middle, p, OUTLET);
//  mgr->connect(m_flow_dir, middle, outlet, p, OUTLET);

  mgr->cellSubdivisionFinished();

  mgr->assignContainer(m_container);
  
  m_container->toVTK(m_geometry_filename);
}
