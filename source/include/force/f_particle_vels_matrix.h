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


#ifndef __F_PARTICLE_VELS_MATRIX_H
#define __F_PARTICLE_VELS_MATRIX_H

#include "f_particle.h"

using namespace std;

class Simulation;

/*!
 * Implementation of a per particle force F acting on the velocities. The force takes a square matrix K as user input and -\a m_factor *K operates on a vector of one user-defined \a point_t U (stored in the tag) of all particles of the user defined colour, i.e. F=-\a m_factor *K.U . This force represents for example a conservative force, if the matrix is a stiffness matrix and the vector U is a displacement, or the force can also represent a friction if the matrix is a friction matrix and the vector U some velocity. NOTE: The force will only work properly if no particles are created or destroyed! It is the user's responsibility to provide a square matrix of appropriate size, which must be in binary format. If addressed by four indices p1, d1, p2, d2 they will be used in exactly this (reverse) order as d2+p2*SPACE_DIMS+d1*SPACE_DIMS*Np+p1*SPACE_DIMS*Np*SPACE_DIMS, where p1, p2 are particle indices, d1, d2 are indexing cartesian coordinates, SPACE_DIMS is the number of dimensions and Np is the number of particles of the user-defined colour.
 */
class FParticleVelsMatrix : public FParticle
{

 private:

  /*!
   * Helper: the size of the matrix
   */
  size_t m_matSize;

  /*!
   * Helper: number of particles (to be assumed constant)
   */
  size_t m_nOfParts;

  protected:

  /*!
   * The matrix provided by the user
   */
  double* m_mat;

  /*!
   * The name of the matrix file provided by the user
   */
  string m_matrixFile;
    
  /*!
   * The name of the input vector the matrix acts on
   */
  string m_invector_name;

  /*!
   * The offset of the input vector the matrix acts on
   */
  size_t m_invec_offset;

  /*!
   * Additional multiplicative factor for the computed force.
   */
  double m_factor;
  
  void init();

 public:
  /*!
   * Constructor
   * @param simulation Pointer to the \a Simulation object
   */
    FParticleVelsMatrix(Simulation *simulation);

  /*!
     * Destructor
   */
    virtual ~FParticleVelsMatrix();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no) {}
#endif

  /*!
     * compute the force
   */
  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP  
    virtual void computeForces(Pairdist* pair, int force_index);
#else
    virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//     virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

  /*!
     * setup this force
   */
    virtual void setup();

  /*!
   * Additional setup
   */
  virtual void setupAfterParticleCreation();

};

#endif
