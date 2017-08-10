/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "particle.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_vector.h"

using namespace std;


#define M_CONTROLLER  ((Controller*) m_parent)
#define M_SIMULATION  ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE  M_SIMULATION->phase()


const Integrator_Register<IntegratorVector> integrator_vector("IntegratorVector");

//---- Constructors/Destructor ----

IntegratorVector::IntegratorVector(Controller *controller): Integrator(controller)
{
  init();
}


IntegratorVector::~IntegratorVector()
{
}



//---- Methods ----

void IntegratorVector::init()
{

  m_properties.setClassName("IntegratorVector");

  m_properties.setDescription(
      "Adds an additional vectorial degree of freedom, to the particles specified. Integration "
      "is performed with a simple Euler scheme."
                             );

  STRINGPC
      (vector, m_vector_name,
       "Full name of the additional vector field, usable as attribute in other modules");

  STRINGPC
      (symbol, m_vector_symbol,
       "Symbol assigned to the additional vector field, usable in algebraic expressions");

  m_vector_name = "vector";
  m_vector_symbol = "sc";
}


void IntegratorVector::setup()
{
  Integrator::setup();

  DataFormat::attribute_t tmpAttr;

  m_vector_offset =
      Particle::s_tag_format[m_colour].addAttribute
      (m_vector_name,
       DataFormat::POINT,
       true,
       m_vector_symbol).offset;

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    tmpAttr =
        Particle::s_tag_format[m_colour].addAttribute
        (STR_FORCE + STR_DELIMITER + m_vector_name + STR_DELIMITER + ObjToString(i),
         DataFormat::POINT,
         true);

    m_force_offset[i] = tmpAttr.offset;
    m_fAttr_index[i] = tmpAttr.index;
  }
  m_dt = M_CONTROLLER->dt();

}

void IntegratorVector::isAboutToStart()
{
  Phase *phase = M_PHASE;

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       for (int j = 0; j < FORCE_HIST_SIZE; ++j)
         i->tag.pointByOffset(((IntegratorVector*) data)->m_force_offset[j]).assign(0);
      );

}

void IntegratorVector::unprotect(size_t index)
{
  Phase *phase = M_PHASE;

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->tag.unprotect(m_fAttr_index[index]);
       if(!index) i->tag.protect(m_fAttr_index[FORCE_HIST_SIZE-1]);
       else i->tag.protect(m_fAttr_index[index-1]);

      );

}



void IntegratorVector::integrateStep1()
{
  Phase *phase = M_PHASE;

  size_t force_index = M_CONTROLLER->forceIndex();

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,

       for(size_t j = 0; j < SPACE_DIMS; ++j)
       {
       i->tag.pointByOffset(((IntegratorVector*) data)->m_vector_offset)[j] +=
           ((IntegratorVector*) data)->m_dt * i->tag.pointByOffset(((IntegratorVector*) data)->m_force_offset[force_index])[j];

      }
      );

}


void IntegratorVector::integrateStep2()
{
}


#ifdef _OPENMP
string IntegratorVector::dofIntegr() {
  return m_vector_name;
}


void IntegratorVector::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->tag.pointByOffset(m_force_offset[force_index])[i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
  }
}

#endif

