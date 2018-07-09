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


#ifndef __F_PARTICLE_VECTOR_H
#define __F_PARTICLE_VECTOR_H

#include "f_particle.h"
#include "function_particle.h"

using namespace std;

class Simulation;

/*!
 * Implementation of a per particle force
 */
class FParticleVector : public FParticle
{
  protected:
  /*!
   * The compiled force expression 
   */
    FunctionParticle m_expression;

  /*!
     * The string for \a m_expression
   */
    string m_exprString;

  /*!
     * The groups this force field should act on. E.g., it is possible to
     * define a force field in the inlet only.
   */
    group_t groups;

      /*!
     * The name of the vector this force acts on
       */
    string m_vector_name;

  /*!
     * The offset of the force on that tensor
   */
    size_t m_force_offset[FORCE_HIST_SIZE];


        
    void init();

  public:
  /*!
   * Constructor
   * @param simulation Pointer to the \a Simulation object
   */
    FParticleVector(Simulation *simulation);

  /*!
     * Destructor
   */
    virtual ~FParticleVector();

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
};

#endif
