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


#include "pc_lattice_bcc.h"

#include "cell.h"
#include "random.h"
#include "simulation.h"
#include "manager_cell.h"


/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorLatticeBCC> particle_creator_lattice_bcc("ParticleCreatorLatticeBCC");


#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())

// bool lattice_defined = false, dd_defined = false;

//---- Constructor/Destructor ----

ParticleCreatorLatticeBCC::ParticleCreatorLatticeBCC(Boundary *boundary): ParticleCreatorLattice(boundary)
{
  init();
}



//---- Methods ----

void ParticleCreatorLatticeBCC::computeDistance()
{
  m_distance = sqrt(3.)/pow(4*m_density, (1./3));
//   m_distance = pow(1/m_density, (1./3));
}

double ParticleCreatorLatticeBCC::spacing(const double& distance) const
{
  return 2*distance/sqrt(3.);
}

void ParticleCreatorLatticeBCC::createParticles()
{    
  //RandomNumberGenerator m_rng;
  point_t box_size, spacing, offset;
  //    cuboid_t box;
  // offset due to wall particles
  //    point_t offset = ((Boundary*) m_parent) -> offset();
  double density;
  //    gsl_rng *rng;

//   size_t seed;

  /* We need the manager to check which group a particle has. */
  ManagerCell *manager = M_PHASE->manager();

  box_size = M_BOUNDARY->boundingBox().size();
  MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "box_size=" << box_size);
  
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

//   if (m_randomize) {
//     seed = 2*getpid() + 1;
//   } else {
//     seed = 1;
//     MSG_INFO("ParticleCreatorLatticeBCC::randomizeVels", "randomize = no --> Numbers won't be random.");
//   }

  //    rng = gsl_rng_alloc(gsl_rng_default);

  //    gsl_rng_set(rng, seed);

  /*    MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "Box size for particles: (" <<
  ((Boundary*) m_parent)->wallBound().x << ", " <<
  ((Boundary*) m_parent)->wallBound().y << ", " <<
  ((Boundary*) m_parent)->wallBound().z << ")");*/

  density = 1;
  for (int i = 0; i < SPACE_DIMS; i++) {/*all this SPECIFIC*/
    
 //     bool dd_is_defined = adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend):dd_defined;
//    spacing[i] = box_size[i]/(m_nlattice_points[i]-1);
    if(m_dd_defined){ /*all this SPECIFIC*/
      spacing[i] = 2*m_distance/(sqrt(3.));
      offset[i] += spacing[i]/4;
      density *= spacing[i];
    }
    else{      /*all this SPECIFIC*/
      // FIXME: the setting of m_distance in this case is a little bit of rubbish and could produce bugs, because m_distance is no array, but spacing is
      int n = m_nlattice_points[i];
      spacing[i] = (double) ((box_size[i] / n ));
      
      m_distance = spacing[i]*sqrt(3.)/2;
      offset[i] += spacing[i]/4;
      density *= spacing[i];
    }
  
  }
    
  MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "spacing = " << spacing);

  /* The density that is generated may differ from the density wished in the configuration
  file because the particle are equally spaced in space and have the same distance from
  the boundaries (of the cube, of course). That means the REAL density can still differ
  from this value */
  
  /*SPECIFIC*/
  density = 2/density;

  MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "density = " << density);

  MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "lattice = " << m_nlattice_points);

  /* All particles sit in a surrounding box of identical size. This is of course assuming the
  Boundary is exactly as big as returnMaxBoxSize proposes. */

  // create the function list for m_properties.unknown()
  // fList will be deleted at the end of this function	
  list<FunctionFixed> fl;
  fList(fl);
  
  // creation of FREE particles
  MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "Creating free particles." << ", spacing = " << spacing << ", offset = " << offset << ", m_nlattice_points = " << m_nlattice_points);
                
  for (int i = 0; i < m_nlattice_points[0]; i++)
    for (int j = 0; j < m_nlattice_points[1]; j++)
      for (int k = 0; k < m_nlattice_points[2]; k++)
        for (int l = 0; l < 2; l++) {
    point_t pos = { { {  /*SPECIFIC*/
      i*spacing[0]+offset[0]+l*spacing[0]/2.,
      j*spacing[1]+offset[1]+l*spacing[1]/2.,
      k*spacing[2]+offset[2]+l*spacing[2]/2.
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
//       else MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "NOT INSIDE: " << p.r << ", i = " << i << ", j = " << j << ", k = " << k << ", l = " << l);
    }
//     else MSG_DEBUG("ParticleCreatorLatticeBCC::createParticles", "NO CELL: " << p.r << ", i = " << i << ", j = " << j << ", k = " << k << ", l = " << l);
    
      }
      scaleVels();

  // next will call the one in ParticleCreatorFree and then ParticleCreatorFreeF::transformVel
      flushParticles();
}

void ParticleCreatorLatticeBCC::init()
{
  m_properties.setClassName("ParticleCreatorLatticeBCC");

  m_properties.setDescription(
      "Generates particles sitting on a body-centred cubic lattice (bcc). This means TWO particles are created per lattice point at the relative positions a*(0.25, 0.25, 0.25) and a*(0.75, 0.75, 0.75), respectively. 'a' is the spacing of the lattice points and will be computed from the given data for 'density', 'distance' and/or 'nLatticeP*' (see below for a description). Note that not all"
      " attributes can be set independently. "
      "For example giving a certain density fixes the distance.\n"
      "Note that the size of the simulation box can be modified by this"
      " ParticleCreator. For example, giving "
      "a density and a number of lattice points fixes the size. If this conflicts with the size "
      "given in the Boundary, the Boundary is being resized to reflect the information given here.\n"
      "The initial conditions of user-defined symbols can be set by taking the symbol as an attribute and defining a mathematical expression for it. For the expression, the same variables are allowed as, e.g., for the attribute 'u'. For non-scalars you have to add the following to the attribute name:\n"
      "\"_x\", \"_y\" or \"_z\" for the respective components of a vector\n"
      "\"_xx\", \"_xy\", \"_xz\", \"_yx\", \"_yy\", \"_yz\", \"_zx\", \"_zy\" and \"_zz\" for the respective components of a tensor."
                             );
    
}


