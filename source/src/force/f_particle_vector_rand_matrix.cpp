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


#include "f_particle_vector_rand_matrix.h"
#include "manager_cell.h"
#include "simulation.h"

using namespace std;

// for SmartEnum
const GenFTypeConcr<FParticleVectorRandMatrix> f_particle_vector_rand_matrix("FParticleVectorRandMatrix");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()
#define M_CONTROLLER M_SIMULATION->controller()


//---- Constructors/Destructor ----

FParticleVectorRandMatrix::FParticleVectorRandMatrix(Simulation *simulation): FParticle(simulation)
{
  init();
} 


FParticleVectorRandMatrix::~FParticleVectorRandMatrix()
{
}


void FParticleVectorRandMatrix::computeForces(Particle* part, int force_index)
{
//   Phase *phase = ((Simulation *) m_parent)->phase();
  
//   size_t nOfParts = phase->returnNofPartC(m_colour);

//   ParticleList& particles = phase->particles(m_colour);

  // slow (constant for this call) particle-row-index
  size_t p1 = part->mySlot;

  /*Check if first particle*/
  if(p1 == 0)
    {
      // Consistency check if it is really the first. In setupAfter... this was also checked, so initially it was OK
      if(((*m_particles)[0]).mySlot != 0)
	throw gError("FParticleVectorRandMatrix::computeForces","Particle with slot 0 is not the first particle anymore. Some ordering or the number of particles have changed. Not supported by this forces! Aborting.");

      /*Create random numbers for this turn*/
      for(size_t p2 = 0; p2 < m_nOfParts; ++p2) {
	for(size_t d2 = 0; d2 < SPACE_DIMS; ++d2) 
	  m_dW[d2+p2*SPACE_DIMS] = m_rng.normal(1);
      }
      //helper
      m_sqrtDt = sqrt(M_CONTROLLER->dt());
    }

  size_t p1Contrib = p1*SPACE_DIMS*m_nOfParts*SPACE_DIMS;

  // the force to be modified
  point_t force = {{{0,0,0}}};

  // outer (slow) cartesian loop
  for(size_t d1 = 0; d1 < SPACE_DIMS; ++d1) {
    size_t d1Contrib = d1*SPACE_DIMS*m_nOfParts;
    // inner (fast) particle loop
    for(size_t p2 = 0; p2 < m_nOfParts; ++p2) {
      // inner (fast) cartesian loop
      size_t p2Contrib = p2*SPACE_DIMS;
      for(size_t d2 = 0; d2 < SPACE_DIMS; ++d2) {
	size_t slot = d2
	  + p2Contrib
	  + d1Contrib
	  + p1Contrib;
	  force[d1] += m_mat[slot]*m_dW[d2+p2Contrib];
      }}}

  // multiply factor and (divide) time-integration correction
  force *= m_factor/m_sqrtDt;

  part->tag.pointByOffset(m_force_offset[force_index]) += force;

// MSG_DEBUG("FParticleVectorRandMatrix::computeForces", "force AFTER = " << part->tag.pointByOffset(m_force_offset[force_index]));

}

#ifndef _OPENMP
void FParticleVectorRandMatrix::computeForces(Pairdist* pair, int force_index)
#else
void FParticleVectorRandMatrix::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FParticleVectorRandMatrix::computeForces", "Fatal error: do not call FParticleVectorRandMatrix::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");  
}


void FParticleVectorRandMatrix::computeForces(int force_index)
{
  throw gError("FParticleVectorRandMatrix::computeForces", "Fatal error: do not call FParticleVectorRandMatrix::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");  
}


void FParticleVectorRandMatrix::init()
{
  m_properties.setClassName("FParticleVectorRandMatrix");

  m_properties.setDescription("Implementation of a per particle force F acting on a vector-symbol. The force takes a square matrix B as user input and factor*B operates on the vector dW, with one random number per particle of the user defined colour, i.e., as SDE, F=factor*K.dW. The corresponding random numbers are drawn from a normal distribution N(mu=0,sigma^2=1). NOTE: The force will only work properly if no particles are created or destroyed! ALSO NOTE that the full matrix is processed, i.e., the algorithm scales with the number of particles squared! It is the user's responsibility to provide a square matrix of appropriate size, which must be in binary format. If the entries are addressed by four indices p1, d1, p2, d2, these indices will be used in exactly this (reverse) order as d2+p2*SPACE_DIMS+d1*SPACE_DIMS*Np+p1*SPACE_DIMS*Np*SPACE_DIMS, where p1, p2 are particle indices, d1, d2 are indexing cartesian coordinates, SPACE_DIMS is the number of dimensions and Np is the number of particles of the user-defined colour. Further note that this Force, as a preparation for Velocity-Verlet time integration, divides its force-contribution by sqrt(dt), with dt the integration-timestep as defined in the Controller.");

  STRINGPC
      (vector, m_vector_name,
       "Name of the vector this force acts on.");

  m_vector_name = "undefined";
   
  STRINGPCINF(nameMatrix, m_matrixFile, "Name of the binary file containing the matrix B.");

  m_matrixFile = "undefined";

  DOUBLEPC(factor, m_factor, -HUGE_VAL, "Additional multiplicative factor for the computed force.");

  m_factor = 1;

  INTPC
    (seed, m_seed, 0,
     "Seed to be used for the random number generator."
     )
    ;
  
  m_seed = RNG_DEFAULT_SEED;

  m_is_pair_force = false;
  m_is_particle_force = true;
  
}

void FParticleVectorRandMatrix::setup()
{

  FParticle::setup();

  m_rng.setSeed(m_seed);

  if(m_vector_name == "undefined")
    throw gError("FParticleVectorRandMatrix::setup", "Attribute 'vector' has value \"undefined\".");
    
  if(m_matrixFile == "undefined")
    throw gError("FParticleVectorRandMatrix::setup", "Attribute 'matrixFile' has value \"undefined\".");
    
    
  DataFormat::attribute_t attr;

  attr = Particle::s_tag_format[m_colour].attrByName(m_vector_name);
  
  
  if(attr.datatype != DataFormat::POINT) 
    throw gError("FParticleVectorRandMatrix::setup", "the symbol " + m_vector_name + 
        " is registered as a non-vector for species " + m_species + "."
                );

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i] =
        Particle::s_tag_format[m_colour].attrByName(string("force_" 
        + m_vector_name + "_" + ObjToString(i))).offset;
  }      

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
    throw gError("FParticleVectorRandMatrix::setup", "Not able to open file " + m_matrixFile + " !");
  }

  M_CONTROLLER->registerForSetupAfterParticleCreation(this);

}

void FParticleVectorRandMatrix::setupAfterParticleCreation()
{
  /*helpers*/
  /*Phase * */  m_phase = ((Simulation *) m_parent)->phase();  
  /*size_t*/ m_nOfParts = m_phase->returnNofPartC(m_colour);
  /*ParticleList* */ m_particles = &(m_phase->particles(m_colour));

  /*reserve memory*/
  m_dW.resize(m_nOfParts*SPACE_DIMS);

  /*Consistency checks*/

  if(((*m_particles)[0]).mySlot != 0)
    throw gError("FParticleVectorRandMatrix::computeForces","Particle with slot 0 is not the first particle. Not supported by this forces! Aborting.");

  if(m_matSize != m_nOfParts*m_nOfParts*SPACE_DIMS*SPACE_DIMS)
    throw gError("FParticleVectorRandMatrix::setupAfterParticleCreation", "Matrix does not seem to have size corresponding to number of particles and dimensions. This is what I have: \n"
		 "Sqrt(matrixSize)/SPACE_DIMS = " + ObjToString(sqrt(m_matSize)/SPACE_DIMS) + "\n"
		 "number of particles = " + ObjToString(m_nOfParts) + "\n"
		 "Aborting."
		 );

}
