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



#include "pc_lattice.h"

#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"


/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorLattice> particle_creator_lattice("ParticleCreatorLattice");


#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())

// bool lattice_defined = false, dd_defined = false;

//---- Constructor/Destructor ----

ParticleCreatorLattice::ParticleCreatorLattice(Boundary *boundary): ParticleCreatorWithRngF(boundary)
{
    init();
}

//---- Methods ----

void ParticleCreatorLattice::computeDistance()
{
  m_distance = pow(1/m_density, (1./3));
}

double ParticleCreatorLattice::spacing(const double& distance) const
{
  return distance;
}

void ParticleCreatorLattice::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend)
{
  m_lattice_defined = false;
  m_dd_defined = false;

  MSG_DEBUG("ParticleCreatorLattice::adjustBoxSize", "size = " << size);
    
  if (m_nlattice_points[0] != -1 &&
      m_nlattice_points[1] != -1 &&
      m_nlattice_points[2] != -1) {
    for (int i = 0; i < SPACE_DIMS; i++)
      if (m_nlattice_points[i] < 1)
        throw gError("ParticleCreatorLattice::adjustBoxSize", "At least one lattice point "
                     "in each direction is needed.");
        
    MSG_DEBUG("ParticleCreatorLattice::adjustBoxSize", "m_nlattice_points = (" 
              << m_nlattice_points[0] << ", " << m_nlattice_points[1] << ", " 
              << m_nlattice_points[2] << ")");
        
    m_lattice_defined = true;
  }
    
  if (m_distance != HUGE_VAL) {
    if (m_distance <= 0)
      throw gError("ParticleCreatorLattice::adjustBoxSize", "Distance must be bigger than zero.");
        
    MSG_DEBUG("ParticleCreatorLattice::adjustBoxSize", "distance = " << m_distance);
        
    m_dd_defined = true;
  }
    
  if (m_density != HUGE_VAL) {
    if (m_dd_defined)
      throw gError("ParticleCreatorLattice::adjustBoxSize", "Please specify only one of: " \
                   "distance, density");
        
    MSG_DEBUG("ParticleCreatorLattice::adjustBoxSize", "density = " << m_density);
      
    /*SPECIFIC*/  
//     m_distance = pow(1/m_density, (1./3));
    
    computeDistance();
    
    MSG_DEBUG("ParticleCreatorLattice::adjustBoxSize", "distance = " << m_distance);        
    m_dd_defined = true;
  }
    
  if (!m_lattice_defined && !m_dd_defined)
    throw gError("ParticleCreatorLattice::adjustBoxSize", "Please specify one of: " \
                 "distance, density, nLatticeP* (with '*' from 'X', 'Y', 'Z')");
					 
    /* If a number of lattice points AND a distance or density AND forceBoundarySize 
    was defined, then the boundary is stretched. */
    if (m_lattice_defined && m_dd_defined && m_force_boundary_size) {
      point_t tmpSize;
      bool mustAbort = false;
      for (int i = 0; i < SPACE_DIMS; i++)
      {
        tmpSize[i] = m_nlattice_points[i] * spacing(m_distance);
        if(tmpSize[i] < 0.999*size[i])
        {
          MSG_INFO("ParticleCreatorLattice::adjustBoxSize", "You are not allowed to"
          " reduce a previously defined box size with settings in this" 
          " ParticleCreatorLattice. This seems to be the case for the " << string(1, 'x'+i) 
          << "-direction:\nold size: " << size[i] << "\nnew size: " << tmpSize[i] 
          << "\nHave to abort soon.");
          mustAbort = true;
        }
      }
      if(mustAbort) throw gError("ParticleCreatorLattice::adjustBoxSize: aborting." 
      " Presumably, you can solve this problem either by defining a smaller box size in"
      " the boundary, or by avoiding two ParticleCreatorLattice to force a box size with" 
      " the second size being smaller than the first.");
      
      size = tmpSize;
      MSG_INFO("ParticleCreatorLattice::adjustBoxSize",
                 "Stretching boundary: size = " << size);
    } 
    else if (!m_lattice_defined) 
    {
        MSG_DEBUG("ParticleCreatorLattice::adjustBoxSize", "Lattice size not defined. "
          "Will be computed, when box size is known.");
        for (int i = 0; i < SPACE_DIMS; i++) {
/*            int n = (int) (size[i] / m_distance + 0.5);
            if (n < 1)
                n = 1;*/
            m_nlattice_points[i] = /*n*/-1;
        }
    } 
}

// new style when inheriting from ParticleCreatorFreeF
void ParticleCreatorLattice::createParticles()
{ 
  //RandomNumberGenerator m_rng;
  point_t box_size, spacing, offset;
  //    cuboid_t box;
  // offset due to wall particles
  //    point_t offset = ((Boundary*) m_parent) -> offset();
  double density;
  //    gsl_rng *rng;

  size_t seed;

  /* We need the manager to check which group a particle has. */
  ManagerCell *manager = M_PHASE->manager();

  box_size = M_BOUNDARY->boundingBox().size();
  MSG_DEBUG("ParticleCreatorLattice::createParticles", "box_size=" << box_size);
  
  if (m_nlattice_points[0] == -1 &&
      m_nlattice_points[1] == -1 &&
      m_nlattice_points[2] == -1) 
  { 
    for (int i = 0; i < SPACE_DIMS; i++) {
      int n = (int) (box_size[i] / m_distance + 0.5);
      if (n < 1)
        n = 1;
      m_nlattice_points[i] = n;
    }
  }

  offset = M_BOUNDARY->boundingBox().corner1;

  initTransform();

  if (m_randomize) {
    seed = 2*getpid() + 1;
  } else {
    seed = 1;
    MSG_INFO("ParticleCreatorLattice::randomizeVels", "randomize = no --> Numbers won't be random.");
  }

  //    rng = gsl_rng_alloc(gsl_rng_default);

  //    gsl_rng_set(rng, seed);

  /*    MSG_DEBUG("ParticleCreatorLattice::createParticles", "Box size for particles: (" <<
  ((Boundary*) m_parent)->wallBound().x << ", " <<
  ((Boundary*) m_parent)->wallBound().y << ", " <<
  ((Boundary*) m_parent)->wallBound().z << ")");*/

  density = 1;
  for (int i = 0; i < SPACE_DIMS; i++) {/*all this SPECIFIC*/
    
 //     bool dd_is_defined = adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend):dd_defined;
//    spacing[i] = box_size[i]/(m_nlattice_points[i]-1);
    if(m_dd_defined){ /*all this SPECIFIC*/
      spacing[i] = m_distance;
      offset[i] += spacing[i]/2;
      density *= spacing[i];
    }
    else{      /*all this SPECIFIC*/
      
      // FIXME: the setting of m_distance in this case is a little bit of rubbish and could produce bugs, because m_distance is no arry, but spacing is
      int n = m_nlattice_points[i];
      m_distance = (double) ((box_size[i] / n ));
      
      spacing[i] = m_distance;
      offset[i] += spacing[i]/2;
      density *= spacing[i];
    }
  
  }
    
  MSG_DEBUG("ParticleCreatorLattice::createParticles", "spacing = " << spacing);

  /* The density that is generated may differ from the density wished in the configuration
  file because the particle are equally spaced in space and have the same distance from
  the boundaries (of the cube, of course). That means the REAL density can still differ
  from this value */
  
  /*SPECIFIC*/
  density = 1/density;

  MSG_DEBUG("ParticleCreatorLattice::createParticles", "density = " << density);

  MSG_DEBUG("ParticleCreatorLattice::createParticles", "lattice = " << m_nlattice_points);

  /* All particles sit in a surrounding box of identical size. This is of course assuming the
  Boundary is exactly as big as returnMaxBoxSize proposes. */

  // create the function list for m_properties.unknown()
  // fList will be deleted at the end of this function	
  list<FunctionFixed> fl;
  fList(fl);
  
  // creation of FREE particles
  MSG_DEBUG("ParticleCreatorLattice::createParticles", "Creating free particles." << ", offset = " << offset);
                
  for (int i = 0; i < m_nlattice_points[0]; i++)
    for (int j = 0; j < m_nlattice_points[1]; j++)
      for (int k = 0; k < m_nlattice_points[2]; k++) {
    point_t pos = { { {  /*SPECIFIC*/
      i*spacing[0]+offset[0],
      j*spacing[1]+offset[1],
      k*spacing[2]+offset[2]
    } } };
    Cell *c;
    Particle p;
	
    p.r = pos;
    p.setColour(m_colour);

    transformPos(p);

    c = manager->findCell(p.r);

    if (c) {
      if (M_BOUNDARY->isInside(p.r)) {
        p.g = c->group();

        for (int w = 0; w < SPACE_DIMS; w++)
          p.v[w] = m_rng.uniform() - 0.5;

        // compute all from m_properties.unknown()
        computeUnknown(fl.begin(), p);
    
        m_particles[p.g].newEntry() = p;
      }
//       else MSG_DEBUG("ParticleCreatorLattice::createParticles", "NOT INSIDE: " << p.r);
    }
      }
      scaleVels();


  // next will call the one in ParticleCreatorFree and then ParticleCreatorFreeF::transformVel
      flushParticles();
}

// old style when inheriting from ParticleCreatorFreePCalc
#if 0
void ParticleCreatorLattice::createParticles()
{    
  //RandomNumberGenerator m_rng;
  point_t box_size, spacing, offset;
  //    cuboid_t box;
  // offset due to wall particles
  //    point_t offset = ((Boundary*) m_parent) -> offset();
  double density;
  //    gsl_rng *rng;

  size_t seed;

  /* We need the manager to check which group a particle has. */
  ManagerCell *manager = M_PHASE->manager();

  box_size = M_BOUNDARY->boundingBox().size();
  MSG_DEBUG("ParticleCreatorLattice::createParticles", "box_size=" << box_size);
  
  if (m_nlattice_points[0] == -1 &&
      m_nlattice_points[1] == -1 &&
      m_nlattice_points[2] == -1) 
  { 
    for (int i = 0; i < SPACE_DIMS; i++) {
      int n = (int) (box_size[i] / m_distance + 0.5);
      if (n < 1)
        n = 1;
      m_nlattice_points[i] = n;
    }
  }

  offset = M_BOUNDARY->boundingBox().corner1;

  initTransform();

  if (m_randomize) {
    seed = 2*getpid() + 1;
  } else {
    seed = 1;
    MSG_INFO("ParticleCreatorLattice::randomizeVels", "randomize = no --> Numbers won't be random.");
  }

  //    rng = gsl_rng_alloc(gsl_rng_default);

  //    gsl_rng_set(rng, seed);

  /*    MSG_DEBUG("ParticleCreatorLattice::createParticles", "Box size for particles: (" <<
  ((Boundary*) m_parent)->wallBound().x << ", " <<
  ((Boundary*) m_parent)->wallBound().y << ", " <<
  ((Boundary*) m_parent)->wallBound().z << ")");*/

  density = 1;
  for (int i = 0; i < SPACE_DIMS; i++) {/*all this SPECIFIC*/
    
 //     bool dd_is_defined = adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend):dd_defined;
//    spacing[i] = box_size[i]/(m_nlattice_points[i]-1);
    if(m_dd_defined){ /*all this SPECIFIC*/
      spacing[i] = m_distance;
      offset[i] += spacing[i]/2;
      density *= spacing[i];
    }
    else{      /*all this SPECIFIC*/
      
      // FIXME: the setting of m_distance in this case is a little bit of rubbish and could produce bugs, because m_distance is no arry, but spacing is
      int n = m_nlattice_points[i];
      m_distance = (double) ((box_size[i] / n ));
      
      spacing[i] = m_distance;
      offset[i] += spacing[i]/2;
      density *= spacing[i];
    }
  
  }
    
  MSG_DEBUG("ParticleCreatorLattice::createParticles", "spacing = " << spacing);

  /* The density that is generated may differ from the density wished in the configuration
  file because the particle are equally spaced in space and have the same distance from
  the boundaries (of the cube, of course). That means the REAL density can still differ
  from this value */
  
  /*SPECIFIC*/
  density = 1/density;

  MSG_DEBUG("ParticleCreatorLattice::createParticles", "density = " << density);

  MSG_DEBUG("ParticleCreatorLattice::createParticles", "lattice = " << m_nlattice_points);

  /* All particles sit in a surrounding box of identical size. This is of course assuming the
  Boundary is exactly as big as returnMaxBoxSize proposes. */

  // creation of FREE particles
  MSG_DEBUG("ParticleCreatorLattice::createParticles", "Creating free particles.");
                
  for (int i = 0; i < m_nlattice_points[0]; i++)
    for (int j = 0; j < m_nlattice_points[1]; j++)
      for (int k = 0; k < m_nlattice_points[2]; k++) {
    point_t pos = { { {  /*SPECIFIC*/
      i*spacing[0]+offset[0],
      j*spacing[1]+offset[1],
      k*spacing[2]+offset[2]
    } } };
    Cell *c;
    Particle p;
	
    p.r = pos;
    p.setColour(m_colour);

    transformPos(p);

    c = manager->findCell(p.r);

    if (c) {
      if (M_BOUNDARY->isInside(p.r)) {
        p.g = c->group();

        for (int w = 0; w < SPACE_DIMS; w++)
          p.v[w] = m_rng.uniform() - 0.5;

        m_particles[p.g].newEntry() = p;
      }
    }
      }
      scaleVels();


  // next will call the one in ParticleCreatorFree
      flushParticles();
}
#endif


void ParticleCreatorLattice::init()
{
  m_properties.setClassName("ParticleCreatorLattice");

  m_properties.setDescription(
    "Generates particles sitting on a simple cubic lattice. This means ONE particle is created per lattice point. Note that not all"
    " attributes can be set independently. "
    "For example giving a certain density fixes the distance.\n"
    "Note that the size of the simulation box can be modified by this"
    " ParticleCreator. For example, giving "
    "a density and a number of lattice points fixes the size. If this conflicts with the size "
      "given in the Boundary, the Boundary is being resized to reflect the information given here.\n" + m_properties.description()
  );
    
  DOUBLEPC
    (kBToverM, m_temperature, 0,
     "k_BT/m = <v^2> for setting the initial thermal velocities of the particles. Note that this "
     "can be overriden by setting vel*.");
    
  for (int i = 0; i < SPACE_DIMS; i++) {
    m_properties.addProperty
      ("nLatticeP" + string(1, 'X'+i), PropertyList::INT, &m_nlattice_points[i],
       new PLCIntGreater(0),
       "Number of lattice points in " + string(1, 'x'+i) + "-"
       "direction.");
    m_nlattice_points[i] = -1;
  }
    
  DOUBLEPC
    (distance, m_distance, 0,
     "Nearest neighbour distance of the particles.");
  DOUBLEPC
    (density, m_density, 0,
     "Number density of the particles.");

  BOOLPC(forceBoundarySize, m_force_boundary_size,
         "Fill the whole Boundary with particles and the density/lattice information "
         "given here. If necessary, the Boundary is going to be resized. You are only " 
         "allowed to INCREASE a previously defined box size.");
    
//   m_randomize = false;
  m_distance = m_density = HUGE_VAL;
  m_temperature = 1;
  m_force_boundary_size = true;
}


void ParticleCreatorLattice::scaleVels()
{
    for (map<int, ParticleList>::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
    //    for (map<int, list<Particle> >::iterator g = m_particles.begin(); g != m_particles.end(); g++) {
        g->second.scaleVels(m_temperature);
    }

}

