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


#include "f_particle_vels_matrix.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<FParticleVelsMatrix> f_particle_vels_matrix("FParticleVelsMatrix");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()
#define M_CONTROLLER M_SIMULATION->controller()

//---- Constructors/Destructor ----

FParticleVelsMatrix::FParticleVelsMatrix(Simulation *simulation): FParticle(simulation)
{
  init();
} 


FParticleVelsMatrix::~FParticleVelsMatrix()
{
}


void FParticleVelsMatrix::computeForces(Particle* part, int force_index)
{
  Phase *phase = ((Simulation *) m_parent)->phase();
  
  size_t nOfParts = phase->returnNofPartC(m_colour);

  ParticleList& particles = phase->particles(m_colour);

  // slow (constant for this call) particle-row-index
  size_t p1 = part->mySlot;

  size_t p1Contrib = p1*SPACE_DIMS*nOfParts*SPACE_DIMS;

  // the force to be modified
  point_t force = {0,0,0};

//   MSG_DEBUG("FParticleVelsMatrix::computeForces", "p=" << p1 << ", force before=" << force);

  // outer (slow) cartesian loop
  for(size_t d1 = 0; d1 < SPACE_DIMS; ++d1) {
    size_t d1Contrib = d1*SPACE_DIMS*nOfParts;
    // inner (fast) particle loop
    for(size_t p2 = 0; p2 < nOfParts; ++p2) {
      point_t& vec = particles[p2].tag.pointByOffset(m_invec_offset);
      // inner (fast) cartesian loop
      for(size_t d2 = 0; d2 < SPACE_DIMS; ++d2) {
	size_t slot = d2
	  + p2*SPACE_DIMS
	  + d1Contrib
	  + p1Contrib;
	  force[d1] += -m_mat[slot]*vec[d2];

// 	  if(m_mat[slot] > 100000000)
// 	    MSG_DEBUG("FParticleVelsMatrix::computeForces", "BUGGY-MAT-ENTRY:p=" << p1 << ", d1=" << d1  << ", p2=" << p2  << ", d2=" << d2 << ", slot=" << slot << ", MAT-ENTRY=" << m_mat[slot]);

// 	  if(vec[d2] > 100000000)
// 	    MSG_DEBUG("FParticleVelsMatrix::computeForces", "BUGGY-VEC-ENTRY:p=" << p1 << ", p2=" << p2  << ", d2=" << d2 << ", VEC-ENTRY=" << vec[d2]);

      }}}

  part->force[force_index] += force*m_factor;


//    MSG_DEBUG("FParticleVelsMatrix::computeForces", "p=" << p1 << ", force after=" << force << ", m_factor=" << m_factor);

}


#ifndef _OPENMP
void FParticleVelsMatrix::computeForces(Pairdist* pair, int force_index)
#else
void FParticleVelsMatrix::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FParticleVelsMatrix::computeForces", "Fatal error: do not call FParticleVelsMatrix::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");  
}


void FParticleVelsMatrix::computeForces(int force_index)
{
  throw gError("FParticleVelsMatrix::computeForces", "Fatal error: do not call FParticleVelsMatrix::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");  
}


void FParticleVelsMatrix::init()
{
  m_properties.setClassName("FParticleVelsMatrix");

  m_properties.setDescription("Implementation of a per particle force F acting on the velocities. The force takes a square matrix K as user input and -factor*K operates on the vector of a user-defined vector-symbol U of all particles of the user defined colour, i.e. F=-factor*K.U . This force represents for example a conservative force, if the matrix is a stiffness matrix and the vector U is a displacement, or the force can also represent a friction if the matrix is a friction matrix and the vector U some velocity. NOTE: The force will only work properly if no particles are created or destroyed! ALSO NOTE that the full matrix is processed, i.e., the algorithm scales with the number of particles squared! It is the user's responsibility to provide a square matrix of appropriate size, which must be in binary format. If the entries are addressed by four indices p1, d1, p2, d2, these indices will be used in exactly this (reverse) order as d2+p2*SPACE_DIMS+d1*SPACE_DIMS*Np+p1*SPACE_DIMS*Np*SPACE_DIMS, where p1, p2 are particle indices, d1, d2 are indexing cartesian coordinates, SPACE_DIMS is the number of dimensions and Np is the number of particles of the user-defined colour.");
  
  STRINGPC
      (inVector, m_invector_name,
       "Name of the input vector the matrix acts on.");

  m_invector_name = "undefined";
 
  STRINGPCINF(nameMatrix, m_matrixFile, "Name of the binary file containing the matrix.");

  m_matrixFile = "undefined";

  DOUBLEPC(factor, m_factor, -HUGE_VAL, "Additional multiplicative factor for the computed force.");

  m_factor = 1;


  m_is_pair_force = false;
  m_is_particle_force = true;
}

void FParticleVelsMatrix::setup()
{


  FParticle::setup();
    
  if(m_invector_name == "undefined")
    throw gError("FParticleVelsMatrix::setup", "Attribute 'inVector' has value \"undefined\".");
    
  if(m_matrixFile == "undefined")
    throw gError("FParticleVelsMatrix::setup", "Attribute 'matrixFile' has value \"undefined\".");
    
    
  DataFormat::attribute_t attr;

  attr = Particle::s_tag_format[m_colour].attrByName(m_invector_name);
    
  if(attr.datatype != DataFormat::POINT) 
    throw gError("FParticleVelsMatrix::setup", "the symbol " + m_invector_name + 
        " is registered as a non-vector for species " + m_species + "."
		 );
 
  m_invec_offset =  Particle::s_tag_format[m_colour].attrByName(m_invector_name).offset;

  // READ matrix
  ifstream file;
  ifstream::pos_type fileSize;

  file.open(m_matrixFile.c_str(), ios::in|ios::binary|ios::ate);
  if (file.is_open()) {
    fileSize = file.tellg();

    m_matSize = (size_t)fileSize/sizeof(double);
    assert(m_matSize%(SPACE_DIMS*SPACE_DIMS) == 0);
    m_mat = new double[m_matSize];

    // go to start and read
    file.seekg (0, ios::beg);
    file.read((char*) m_mat, fileSize);
    file.close();
  }
  else {
    throw gError("FParticleVelsMatrix::setup", "Not able to open file " + m_matrixFile + " !");
  }

  M_CONTROLLER->registerForSetupAfterParticleCreation(this);

}

void FParticleVelsMatrix::setupAfterParticleCreation()
{
  /*helpers*/
  Phase*  m_phase = ((Simulation *) m_parent)->phase();  
  size_t  m_nOfParts = m_phase->returnNofPartC(m_colour);
  ParticleList& m_particles = m_phase->particles(m_colour);

  /*Consistency checks*/

  if(m_particles[0].mySlot != 0)
    throw gError("FParticleVectorRandMatrix::computeForces","Particle with slot 0 is not the first particle. Not supported by this forces! Aborting.");

  if(m_matSize != m_nOfParts*m_nOfParts*SPACE_DIMS*SPACE_DIMS)
    throw gError("FParticleVectorRandMatrix::setupAfterParticleCreation", "Matrix does not seem to have size corresponding to number of particles and dimensions. This is what I have: \n"
		 "Sqrt(matrixSize)/SPACE_DIMS = " + ObjToString(sqrt(m_matSize)/SPACE_DIMS) + "\n"
		 "number of particles = " + ObjToString(m_nOfParts) + "\n"
		 "Aborting."
		 );

}
