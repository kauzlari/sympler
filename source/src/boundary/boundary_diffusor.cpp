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



#include "boundary_diffusor.h"

#include "manager_cell.h"
#include "wall_triangle.h"
#include "simulation.h"
#include "pc_inlet.h"


/* Register this Boundary with the factory. */
const Boundary_Register<BoundaryDiffusor> boundary_diffusor("BoundaryDiffusor");


#define FRAME_OFF 0

#define FLOW_DIR 2
#define PERP_DIR1 0
#define PERP_DIR2 1

#define M_PHASE  ((Phase*) m_parent)

//---- Constructors/Destructor ----

BoundaryDiffusor::BoundaryDiffusor(Phase *phase): BoundaryWithInlet(phase)
{
	init();
}


BoundaryDiffusor::~BoundaryDiffusor()
{
}



//---- Methods ----

void BoundaryDiffusor::init()
{
  m_properties.setClassName("BoundaryWithInlet");
  m_properties.setName("BoundaryDiffusor");

  m_properties.setDescription
    ("Simple diffusor type geometry. Cuboid inlet which expands to a cuboid reservoir.\n\n"
      

      "                                      .--------------------.\n"
      "                                     .                    .|\n"
      "           .------------.-----------.                    . |\n"
      "          .            .           .                    .  |\n"
      " IWidth  .            .            --------------------    |\n"
      "        .            .            |                    |   | DDepth\n"
      "        --------------------------|                    |   |\n"
      "       |            :             |                    |   |\n"
      "       |  inlet (I) :             |     diffusor (D)   |   |\n"
      "IDepth |            :             |                    |   .\n"
      "       |            :             |                    |  .\n"
      "        --------------------------|                    | . DWidth\n"
      "           ILength  :   ILength   |                    |.\n"
      "y ^                                --------------------\n"
      "  |  x/                                  DLength\n"
      "  |  /\n"
      "  | /\n"
      "   ----->\n"
      "       z\n"
     );

  BOOLPC
      (inlet, m_inlet,
       "Should we have a pressure inlet?");
  BOOLPC
      (periodic, m_periodic,
       "Should the y-direction be periodic?");
  DOUBLEPC
      (inletWidth, m_inlet_width, 0,
       "Width of the inlet (see diagram above).");
  DOUBLEPC
      (inletDepth, m_inlet_depth, 0,
       "Depth of the inlet (see diagram above).");
/*  DOUBLEPC
    (inletHeight, m_inlet_height, 0,
     "Height of the inlet (see diagram above).");*/
  DOUBLEPC
      (diffusorWidth, m_diffusor_width, 0,
       "Width of the diffusor (see diagram above).");
  DOUBLEPC
      (diffusorDepth, m_diffusor_depth, 0,
       "Depth of the diffusor (see diagram above).");
  DOUBLEPC
    (diffusorLength, m_diffusor_length, 0,
     "Length of the diffusor in z-direction (see diagram above).");

  m_periodic = true;
  m_inlet_width = 1;
  m_inlet_depth = 1;
//   m_inlet_height = 2;
  m_diffusor_width = 2;
  m_diffusor_depth = 1;
  m_diffusor_length = 2;
}

#define ADDRECT(a, b, c, d)  {						\
    /*MSG_DEBUG("BoundaryDiffusor::read", "Adding: " << a << ", " << b << ", " << c << ", " << d);*/ \
    m_container->addWall(new WallTriangle(m_container, m_container->reflector(), a, b, c)); \
    m_container->addWall(new WallTriangle(m_container, m_container->reflector(), c, d, a)); } while(0)
	
	// Was a, b, c | b, c, d
	
#define ADDVERTEX(v, a, b, c) {			\
    point_t p;					\
    p.x = a+FRAME_OFF;				\
    p.y = b+FRAME_OFF;				\
    p.z = c+FRAME_OFF;				\
    v = m_container->addVertex(p); } while(0)

#define ADDRECTTOBASE(a, b, c, d)  {						\
   m_inlet_base_surface.addWall(new WallTriangle(m_container, m_inlet_base_surface.reflector(), a, b, c)); \
    m_inlet_base_surface.addWall(new WallTriangle(m_container, m_inlet_base_surface.reflector(), c, d, a)); } while(0)
	
	// Was a, b, c | b, c, d
	
#define ADDVERTEXTOBASE(v, a, b, c) {			\
    point_t p;					\
    p.x = a+FRAME_OFF;				\
    p.y = b+FRAME_OFF;				\
    p.z = c+FRAME_OFF;				\
    v = m_inlet_base_surface.addVertex(p); } while(0)


void BoundaryDiffusor::setup()
{
  int v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;
  int bv0, bv1, bv2, bv3;
  double help, helpDepth;

  BoundaryArbitrary::setup();

  m_container = new WallContainer(findReflector("default"));

  help = m_diffusor_width/2 - m_inlet_width/2;
  if(help < 0) throw gError("BoundaryDiffusor::setup", "Currently, inletWidth > diffusorWidth is not allowed!");
  
  helpDepth = m_diffusor_depth/2 - m_inlet_depth/2;
  if(helpDepth < 0) throw gError("BoundaryDiffusor::setup", "Currently, inletDepth > diffusorDepth is not allowed!");

  ADDVERTEX(v0, help, helpDepth, 0);
  ADDVERTEX(v1, help, helpDepth+m_inlet_depth, 0);
  ADDVERTEX(v2, help+m_inlet_width, helpDepth+m_inlet_depth, 0);
  ADDVERTEX(v3, help+m_inlet_width, helpDepth, 0);

  ADDVERTEXTOBASE(bv0, help, helpDepth, 0);
  ADDVERTEXTOBASE(bv1, help, helpDepth+m_inlet_depth, 0);
  ADDVERTEXTOBASE(bv2, help+m_inlet_width, helpDepth+m_inlet_depth, 0);
  ADDVERTEXTOBASE(bv3, help+m_inlet_width, helpDepth, 0);

  ADDVERTEX(v4, help, helpDepth, 2*m_inlet_length);
  ADDVERTEX(v5, help, helpDepth+m_inlet_depth, 2*m_inlet_length);
  ADDVERTEX(v6, help+m_inlet_width, helpDepth+m_inlet_depth, 2*m_inlet_length);
  ADDVERTEX(v7, help+m_inlet_width, helpDepth, 2*m_inlet_length);

  ADDVERTEX(v8, 0, 0, 2*m_inlet_length);
  ADDVERTEX(v9, 0, m_diffusor_depth, 2*m_inlet_length);
  ADDVERTEX(v10, m_diffusor_width, m_diffusor_depth, 2*m_inlet_length);
  ADDVERTEX(v11, m_diffusor_width, 0, 2*m_inlet_length);

  ADDVERTEX(v12, 0, 0, 2*m_inlet_length+m_diffusor_length);
  ADDVERTEX(v13, 0, m_diffusor_depth, 2*m_inlet_length+m_diffusor_length);
  ADDVERTEX(v14, m_diffusor_width, m_diffusor_depth, 2*m_inlet_length+m_diffusor_length);
  ADDVERTEX(v15, m_diffusor_width, 0, 2*m_inlet_length+m_diffusor_length);

  ADDRECT(v3, v2, v1, v0);
  ADDRECTTOBASE(bv3, bv2, bv1, bv0);
  
  if(!m_periodic)
  {
    ADDRECT(v0, v4, v7, v3);
    ADDRECT(v2, v6, v5, v1);
  }
  // consistency check for periodic case
  else
  {
    if(helpDepth != 0.)
      throw gError("BoundaryDiffusor::setup", "Periodic boundary conditions only possible if diffusorDepth = inletDepth.");
  }
  
  ADDRECT(v3, v7, v6, v2);
  ADDRECT(v1, v5, v4, v0);

  if(!m_periodic)
  {
    ADDRECT( v8, v12, v15, v11);
    ADDRECT(v10, v14, v13,  v9);
  }

    
  ADDRECT(v11, v15, v14, v10);
  ADDRECT( v9, v13, v12,  v8);

  ADDRECT(v12, v13, v14, v15);

  // if the helpers are zero, the corresponding quadrilaterals would have zero area
  if(helpDepth > 0)
  {
    ADDRECT( v8, v11,  v7,  v4);
    ADDRECT(v10,  v9,  v5,  v6);
  }
  if(help > 0)
  {
    ADDRECT(v11, v10,  v6,  v7);
    ADDRECT( v9,  v8,  v4,  v5);
  }
  
  m_container->updateBoundingBox();
	
  MSG_DEBUG("BoundaryDiffusor::setup", "box: " << m_container->boundingBox().corner1 << " - " << m_container->boundingBox().corner2);
}



//---- Cell subdivision ----

void BoundaryDiffusor::setup(Simulation* sim, ManagerCell *mgr)
{
  point_t old_size;
  old_size = m_container->boundingBox().corner2 - m_container->boundingBox().corner1;
	
  /*	m_size = m_container->boundingBox().corner2;
    old_size = m_size;

    MSG_DEBUG("BoundaryDiffusor::setup", "old_size = " << old_size);

    // m_particle_creator->adjustBoxSize(m_size);

    bool_point_t frontFrame, endFrame;
    m_particle_creator->adjustBoxSize(m_size, frontFrame, endFrame);


    m_volume = m_size.x * m_size.y * m_size.z; */

  // periodicity
  bool_point_t periodicity;

  periodicity.x = periodicity.z = false;
  
  periodicity.y = m_periodic;
  m_container->setPeriodicityFront(periodicity);
  m_container->setPeriodicityBack(periodicity);
  
  // now the parent class setups can let the ParticleCreators adjust the box size and the frame 
  BoundaryWithInlet::setup(sim, mgr);
	
  if (!(old_size == m_proposedSize)) {
    throw gError("BoundaryDiffusor::setup", "Cannot resize a diffusor.");
  }
	

  /* Cell subdivision */
  point_t c1, c2;
  region_t *inlet, *diffusor;
//   bool_point_t p;
	
  MSG_DEBUG("BoundaryDiffusor::setup", "Initalizing cells.");

  c1 = m_container->boundingBox().corner1;
  c2[1] = m_diffusor_depth;

  if((m_frontFrame[1] || m_endFrame[1]) && m_periodic)
    throw gError("BoundaryDiffusor::setup", "cannot have wall particles at a periodic face");
  if(m_endFrame[1]) c2[1] += M_PHASE -> pairCreator() -> interactionCutoff();
//   if(m_endFrame[1]) c2[1] += sim -> maxCutoff;
  
  c2[0] = m_diffusor_width;

  if(m_endFrame[0]) c2[0] += M_PHASE -> pairCreator() -> interactionCutoff();
//   if(m_endFrame[0]) c2[0] += sim -> maxCutoff;
 
  if(m_inlet)
  {
    c2[2] = m_inlet_length;
    ParticleCreatorInlet *pc = NULL;

    for(vector<ParticleCreator*>::iterator pcIter = m_pcList.begin();
        pcIter != m_pcList.end(); ++pcIter)
    {
      if(typeid(**pcIter) == typeid(ParticleCreatorInlet))
      {
        if(pc) 
          throw gError("BoundarySTL::setup", "Please define only one ParticleCreatorInlet");
        pc = (ParticleCreatorInlet*) (*pcIter);
        MSG_DEBUG("BoundaryDiffusor::setup", "Adding inlet.");
      }
    }

    if(!pc) throw gError("BoundarySTL::setup", "Please define a ParticleCreatorInlet.");
    		   
    // will be needed by the ParticleCreatorInlet
    m_inlet_normal.x = m_inlet_normal.y = 0;
    m_inlet_normal.z = -1;
    
    inlet = mgr->cellSubdivide(sim->maxCutoff, c1, c2, periodicity, 1, pc, 2, 1, P_CREATE);
    
    pc->setRegion(inlet);
      
    c1[2] = c2[2];
    c2[2] = 2*m_inlet_length+m_diffusor_length;
    if(m_endFrame[2]) c2[2] += M_PHASE -> pairCreator() -> interactionCutoff();
//     if(m_endFrame[2]) c2[2] += sim -> maxCutoff;
    
    diffusor = mgr->cellSubdivide(sim->maxCutoff, c1, c2, periodicity, 0, pc, 2, -1, P_DELETE);

    MSG_DEBUG("BoundaryDiffusor::setup", "Establishing connections between inlet and diffusor.");

            
    /* Connect inlet to diffusor */
    
    ManagerCell::connect(2, inlet, diffusor, periodicity, NEIGHBOR);
    
  }
  else
  {
    c2[2] = 2*m_inlet_length+m_diffusor_length;
    if(m_endFrame[2]) c2[2] += M_PHASE -> pairCreator() -> interactionCutoff();
//     if(m_endFrame[2]) c2[2] += sim -> maxCutoff;
    diffusor = mgr->cellSubdivide(sim->maxCutoff, c1, c2, periodicity, 0);
  }
                     
  mgr->cellSubdivisionFinished();

  mgr->assignContainer(m_container);
  
  m_container->toVTK(m_geometry_filename);
}

