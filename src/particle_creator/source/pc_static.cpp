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




#include "pc_static.h"

#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"


/*!
 * Register this ParticleCreator with the factory. 
 */
const ParticleCreator_Register<ParticleCreatorStatic> particle_creator_static("ParticleCreatorStatic");


#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())

// bool static_defined = false, ddstatic_defined = false;

/*!
 *---- Constructor/Destructor ----
 */

ParticleCreatorStatic::ParticleCreatorStatic(Boundary *boundary): ParticleCreatorWithRngPCalc(boundary)
{
  init();
}


void ParticleCreatorStatic::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend)
{
  m_static_defined = false;
  m_ddstatic_defined = false;

  MSG_DEBUG("ParticleCreatorStatic::adjustBoxSize", "size = " << m_static_box.size());
  
  /*
   * You can only define nStaticPoints, density or distance
   */
    
  if (m_nstatic_points[0] != -1 &&
      m_nstatic_points[1] != -1 &&
      m_nstatic_points[2] != -1) {
      
    if(m_ddstatic_defined)
      throw gError("ParticleCreatorStatic::adjustBoxSize", "You can only choose one of density,"
                     "distance, nStaticPX/Y/Z.");
      
    for (int i = 0; i < SPACE_DIMS; i++)
      if (m_nstatic_points[i] < 1)
        throw gError("ParticleCreatorStatic::adjustBoxSize", "At least one static point "
                     "in each direction is needed.");
       
    MSG_DEBUG("ParticleCreatorStatic::adjustBoxSize", "m_nstatic_points = (" 
              << m_nstatic_points[0] << ", " << m_nstatic_points[1] << ", " 
              << m_nstatic_points[2] << ")");
        
    m_static_defined = true;
  }
    
  if (m_distance != HUGE_VAL) {
    
    if (m_static_defined || m_ddstatic_defined)
      throw gError("ParticleCreatorStatic::adjustBoxSize", "You cannot define distance, density and nStaticP at once.");
    
      if (m_distance <= 0)
        throw gError("ParticleCreatorStatic::adjustBoxSize", "Distance must be bigger than zero.");
      
    MSG_DEBUG("ParticleCreatorStatic::adjustBoxSize", "distance = " << m_distance);
        
    m_ddstatic_defined = true; 
    
  }

    
  if (m_density != HUGE_VAL) {
    if (m_ddstatic_defined || m_static_defined)
      throw gError("ParticleCreatorStatic::adjustBoxSize", "Please specify only one of: " \
                   "distance, density, nStaticP");
        
    MSG_DEBUG("ParticleCreatorStatic::adjustBoxSize", "density = " << m_density);
        
    m_distance = pow(1/m_density, (1./3));
    MSG_DEBUG("ParticleCreatorStatic::adjustBoxSize", "distance = " << m_distance);        
    m_ddstatic_defined = true;
  }
    
  if (!m_static_defined && !m_ddstatic_defined)
    throw gError("ParticleCreatorStatic::adjustBoxSize", "Please specify one of: " \
                 "distance, density, nStaticP* (with '*' from 'X', 'Y', 'Z')");
 
}

/*! 
 *createParticles: Almost like ParticleCreatorLattice with some modifications (see PCLattice)
 */
 
void ParticleCreatorStatic::createParticles()
{ 
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "corner1=" << m_static_box.corner1);
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "corner2=" << m_static_box.corner2);
      point_t box_size, spacing;
      //   size_t seed;
      double density;
      ManagerCell *manager;
      m_offset = m_static_box.corner1;
      box_size = m_static_box.size();
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "m_offset=" << m_offset);
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "box_size=" << box_size);
  
      /*
       *When the staticBox is bigger than the boundingBox or outside of it ,there'll come out an error. This checks these cases.
       */
   
      point_t c1 = M_BOUNDARY->boundingBox().corner1;
      point_t c2 = M_BOUNDARY->boundingBox().corner2;
      point_t c11 = m_static_box.corner1;
      point_t c22 = m_static_box.corner2;
      point_t boxcenter = (m_static_box.corner1 + m_static_box.corner2)/2;
      point_t boxradii = (m_static_box.corner2 - m_static_box.corner1)/2;
      for(int i = 0; i < SPACE_DIMS; i++)
	boxradii[i] = abs(boxradii[i]);
           
      manager = M_PHASE->manager();
      
      if (m_nstatic_points[0] == -1 &&
	  m_nstatic_points[1] == -1 &&
	  m_nstatic_points[2] == -1) 
	{ 
	  for (int i = 0; i < SPACE_DIMS; i++) {
	    int n = (int) (box_size[i] / m_distance + 0.5);
	    if (n < 1)
	      n = 1;
	    m_nstatic_points[i] = n;
	  }
	}
      
      initTransform();
            
      density = 1;
      for (int i = 0; i < SPACE_DIMS; i++) {

	if(m_ddstatic_defined){
	  spacing[i] = m_distance;
	  m_offset[i] += spacing[i]/2;
	  density *= spacing[i];
	}
	else {
	  int n = m_nstatic_points[i];
	  m_distance = (double) ((box_size[i] / n ));
	  
	  spacing[i] = m_distance;
	  m_offset[i] += spacing[i]/2;
	  density *= spacing[i];
	}
      }
      
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "spacing = " << spacing);
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "offset = " << m_offset);
  
  density = 1/density;
  
  MSG_DEBUG("ParticleCreatorStatic::createParticles", "density = " << density);
  
  MSG_DEBUG("ParticleCreatorStatic::createParticles", "Static corner1 = " << m_static_box.corner1);
  
  MSG_DEBUG("ParticleCreatorStatic::createParticles", "Static corner2 = " << m_static_box.corner2);

  MSG_DEBUG("ParticleCreatorStatic::createParticles", "# static points = " << m_nstatic_points);

  MSG_DEBUG("ParticleCreatorStatic::createParticles", "Creating free particles.");
  
  if(m_timeOffset == 0. || !m_more)
    {
      // the next one is needed if we create more particles later on (m_more = true)
      m_nextTimeForMore = m_dtForMore; 

      /* 
       * Creates the particles according to the defined nStaticP's
       */
      
      for (int i = 0; i < m_nstatic_points[0]; i++)
	for (int j = 0; j < m_nstatic_points[1]; j++)
	  for (int k = 0; k < m_nstatic_points[2]; k++) {
	    point_t pos = { { {
		  i*spacing[0]+m_offset[0],
		  j*spacing[1]+m_offset[1],
		  k*spacing[2]+m_offset[2]
		} } };
	    Cell *c;
	    Particle p;
	    
	    p.r = pos;
	    p.setColour(m_colour);
	    
	    transformPos(p);
	    
	    c = manager->findCell(p.r);
	    
	    if (c) {
	      if (!m_checkInside || M_BOUNDARY->isInside(p.r)) {
		if(!m_ellipsoid || ((p.r-boxcenter).listDivide(boxradii)).absSquare() <= 1) {
		  p.g = c->group();
		  
		  for (int w = 0; w < SPACE_DIMS; w++)
		    p.v[w] = m_rng.uniform() - 0.5;
		  
		  m_particles[p.g].newEntry() = p;
		}
// 		else MSG_DEBUG ("ParticleCreatorStatic::createParticles", "elipsoid check failed: " << p.r);
	      }
// 	      else {MSG_DEBUG ("ParticleCreatorStatic::createParticles", "not inside: " << p.r);}
	      
	    }
// 	     	else {MSG_DEBUG ("ParticleCreatorStatic::createParticles", "no cell, r=" << p.r);}
	    
	    //         cout.precision(21);
	    //         MSG_DEBUG("ParticleCreatorStatic::createParticle", "pos = " << p.r);
	  }
      
      scaleVels();
      
      flushParticles();
    }
  else
    {
      m_nextTimeForMore = m_timeOffset; // so nothing created now
      MSG_DEBUG("ParticleCreatorStatic::createParticles", "NO INITIAL CREATION, m_nextTimeForMore = " << m_nextTimeForMore);
    }
}

void ParticleCreatorStatic::createMoreParticles()
{
  //  MSG_DEBUG("ParticleCreatorStatic::createMoreParticles", "START");
  if(m_more)
  {
    //     MSG_DEBUG("ParticleCreatorStatic::createMoreParticles", "time-check: simtime = " << ((Simulation*)(M_PHASE->parent()))->controller()->time() << ", nextTime = " << m_nextTimeForMore);
    if
    (
      m_nextTimeForMore <=
      (
        (Simulation*)
          (M_PHASE->parent())
      )
        ->controller()->time()
    )
    {
      m_nextTimeForMore += m_dtForMore;
      MSG_DEBUG("ParticleCreatorStatic::createMoreParticles", "creating new particles. Next time for creation: " << m_nextTimeForMore);
      point_t box_size, spacing;
      double density;
      ManagerCell *manager;
      m_offset = m_static_box.corner1;
      box_size = m_static_box.size();
      point_t boxcenter = (m_static_box.corner1 + m_static_box.corner2)/2;
      point_t boxradii = (m_static_box.corner2 - m_static_box.corner1)/2;
      for(int i = 0; i < SPACE_DIMS; i++)
    	  boxradii[i] = abs(boxradii[i]);

      manager = M_PHASE->manager();
    
      // next is done once in createParticles
    //   if (m_nstatic_points[0] == -1 &&
    //       m_nstatic_points[1] == -1 &&
    //       m_nstatic_points[2] == -1) 
    //   { 
    //     for (int i = 0; i < SPACE_DIMS; i++) {
    //       int n = (int) (box_size[i] / m_distance + 0.5);
    //       if (n < 1)
    //         n = 1;
    //       m_nstatic_points[i] = n;
    //     }
    //   }
    
    //   initTransform();
    
    //   if (m_randomize) {
    //     seed = 2*getpid() + 1;
    //   } else {
    //     seed = 1;
    //     MSG_INFO("ParticleCreatorStatic::randomizeVels", "randomize = no --> Numbers won't be random.");
    //   }
    
      density = 1;
      for (int i = 0; i < SPACE_DIMS; ++i) {
      
        if(m_ddstatic_defined){
          spacing[i] = m_distance;
          m_offset[i] += spacing[i]/2;
          density *= spacing[i];
        }
        else {
          int n = m_nstatic_points[i];
          m_distance = (double) ((box_size[i] / n ));
          
          spacing[i] = m_distance;
          m_offset[i] += spacing[i]/2;
          density *= spacing[i];
        }
      }
        
      /*! 
      * Creates the particles according to the defined nStaticP's
      */
    
      for (int i = 0; i < m_nstatic_points[0]; i++)
        for (int j = 0; j < m_nstatic_points[1]; j++)
          for (int k = 0; k < m_nstatic_points[2]; k++) {
        point_t pos = { { {
          i*spacing[0]+m_offset[0],
          j*spacing[1]+m_offset[1],
          k*spacing[2]+m_offset[2]
        } } };
        Cell *c;
        Particle p;
        Particle *new_p;
      
        p.r = pos;
        p.setColour(m_colour);
    
        transformPos(p);
    
        c = manager->findCell(p.r);
    
        if (c) {
          if (!m_checkInside || M_BOUNDARY->isInside(p.r)) {
        	  if(!m_ellipsoid || ((p.r-boxcenter).listDivide(boxradii)).absSquare() <= 1) {
                p.g = c->group();
    
    /*        for (int w = 0; w < SPACE_DIMS; w++)
              p.v[w] = m_rng.uniform() - 0.5;
    
            m_particles[p.g].newEntry() = p;*/
          
    //       else {MSG_DEBUG ("ParticleCreatorStatic::createParticles", "not inside");}
    
              p.dt = 0;
      
      //       scaleVels();
      
      //       flushParticles();
        
              for (int w = 0; w < SPACE_DIMS; w++)
                p.v[w] = m_rng.normal(m_temperature_sqrt);
//             MSG_DEBUG ("ParticleCreatorStatic::createMoreParticles", "BEFORE transfromVel: p.v = " << p.v);      
              transformVel(p);
//             MSG_DEBUG ("ParticleCreatorStatic::createMoreParticles", "AFTER transfromVel: p.v = " << p.v);      
        
              new_p = M_PHASE->addParticle(p);
//               MSG_DEBUG("ParticleCreatorStatic::createParticle", "newp: " << new_p->mySlot/*r << endl << new_p->v*/);
              c->injectFree(m_colour, new_p);
            }
          }
        }
        else 
        {
          throw gError("ParticleCreatorStatic::createParticles", "no cell for particle with r = (" + ObjToString(p.r.x) + ", " + ObjToString(p.r.y) + ", " + ObjToString(p.r.z) + ").");
        }
      }
    }
  }
}

/*!
 * PC_STATIC init
*/

void ParticleCreatorStatic::init()
{
  m_properties.setClassName("ParticleCreatorStatic");
  
  m_properties.setDescription
    ("Creates particles in a simple cubic lattice inside a user-defined cuboid box that not necessarily coincides with the simulation boundaries. As an option, the particles may be constructed inside an ellipsoid.");
  
  BOOLPC
    (more, m_more,
     "Should more particles be created DURING the simulation?");
  
  BOOLPC
    (checkIfInside, m_checkInside,
     "Should the particle creator explicitly check whether the particles that should be created are really inside and only create them if this is true?");
  
  DOUBLEPC
    (temperature, m_temperature, 0,
     "Initial temperature of the particles. Note that this "
     "can be overriden by setting vel*.");
  
  
  DOUBLEPC
    (dtForMore, m_dtForMore, 0,
     "Time between creation of more particles. Will be ignored if 'more = \"no\"'.");
  
  DOUBLEPC
    (timeOffset, m_timeOffset, 0,
     "Time offset for the particle creation. The first particles are created at the time t = timeOffset. Will be ignored if 'more = \"no\"'.");
  
  
  for (int i = 0; i < SPACE_DIMS; i++) {
    m_properties.addProperty
      ("nStaticP" + string(1, 'X'+i), PropertyList::INT, &m_nstatic_points[i],
       new PLCIntGreater(0),
       "Number of particles created in " + string(1, 'x'+i) + "-"
       "direction.");
    m_nstatic_points[i] = -1;
  }
    
  m_properties.addProperty
    ("corner1", PropertyList::POINT, &m_static_box.corner1,
     NULL,
     "Sets the lower limits of the cuboid.");
  
  m_properties.addProperty
    ("corner2", PropertyList::POINT, &m_static_box.corner2,
     NULL,
     "Sets the upper limits of the cuboid.");
  
  DOUBLEPC
    (distance, m_distance, 0,
     "Spacing of particles.");
  DOUBLEPC
    (density, m_density, 0,
     "Density of the particles.");
  
  BOOLPC(forceBoundarySize, m_force_boundary_size,
         "Fill the whole Boundary with particles and the density/lattice information "
         "given here. If necessary, the Boundary is going to be resized. You are only " 
         "allowed to INCREASE a previously defined box size.");
  
  BOOLPC(ellipsoid, m_ellipsoid,
         "The particles form the inscribed ellipsoid of the box given by the corners. "
         "Useful for creating droplets.");
  
  /* Default values */
    
  m_distance = m_density = HUGE_VAL;
  m_temperature = 1;
  m_force_boundary_size = true;
  
  m_static_box.corner1[0] = 0;
  m_static_box.corner1[1] = 0;
  m_static_box.corner1[2] = 0;
  
  m_static_box.corner2[0] = 0;
  m_static_box.corner2[1] = 0;
  m_static_box.corner2[2] = 0;

  for (int i = 0; i < SPACE_DIMS; i++) {
    if ((m_static_box.corner1[i] == m_static_box.corner2[i]) == 0) 
      {
        throw gError("ParticleCreatorStatic::init", "The box-size that you've defined is 0. Please correct this.");
      }
  }

  m_more = false;  
  m_checkInside = true;
  m_ellipsoid = false;

//   m_dtForMore.addVariable("t");
//   m_dtForMore.setExpression("1");
  m_dtForMore = 1;
  m_timeOffset = 0;
}

void ParticleCreatorStatic::scaleVels()
{
  for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
    g->second.scaleVels(m_temperature);
  }
}
