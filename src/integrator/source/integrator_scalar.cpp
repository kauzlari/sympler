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
#include "integrator_scalar.h"

using namespace std;


#define M_CONTROLLER  ((Controller*) m_parent)
#define M_SIMULATION  ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE  M_SIMULATION->phase()


const Integrator_Register<IntegratorScalar> integrator_scalar("IntegratorScalar");

//---- Constructors/Destructor ----

IntegratorScalar::IntegratorScalar(Controller *controller): Integrator(controller)
{
  init();
}


IntegratorScalar::~IntegratorScalar()
{
}


//---- Methods ----

void IntegratorScalar::init()
{

  m_properties.setClassName("IntegratorScalar");

  m_properties.setDescription(
    "Adds an additional degree of freedom, to the particles specified. Integration "
      "is performed with a simple Euler scheme."
  );

  STRINGPC
    (scalar, m_scalar_name,
     "Full name of the additional scalar field, usable as attribute in other modules");

  STRINGPC
      (symbol, m_scalar_symbol,
       "Symbol assigned to the additional scalar field, usable in algebraic expressions");

  m_scalar_name = "scalar";
  m_scalar_symbol = "sc";
}


void IntegratorScalar::setup()
{
  Integrator::setup();

  DataFormat::attribute_t tmpAttr;

  m_scalar_offset =
    Particle::s_tag_format[m_colour].addAttribute
      (m_scalar_name,
       DataFormat::DOUBLE,
       true,
       m_scalar_symbol).offset;

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    tmpAttr =
        Particle::s_tag_format[m_colour].addAttribute
        (STR_FORCE + STR_DELIMITER + m_scalar_name + STR_DELIMITER + ObjToString(i),
         DataFormat::DOUBLE,
         true);

    m_force_offset[i] = tmpAttr.offset;
    m_fAttr_index[i] = tmpAttr.index;
  }
  m_dt = M_CONTROLLER->dt();

}

void IntegratorScalar::isAboutToStart()
{
  Phase *phase = M_PHASE;

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       for (int j = 0; j < FORCE_HIST_SIZE; ++j)
         i->tag.doubleByOffset(((IntegratorScalar*) data)->m_force_offset[j]) = 0;
      );

}

void IntegratorScalar::unprotect(size_t index)
{
  Phase *phase = M_PHASE;

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->tag.unprotect(m_fAttr_index[index]);
       if(!index) i->tag.protect(m_fAttr_index[FORCE_HIST_SIZE-1]);
       else i->tag.protect(m_fAttr_index[index-1]);

      );

}



void IntegratorScalar::integrateStep1()
{
  Phase *phase = M_PHASE;

  size_t force_index = M_CONTROLLER->forceIndex();

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,

       i->tag.doubleByOffset(((IntegratorScalar*) data)->m_scalar_offset) +=
           ((IntegratorScalar*) data)->m_dt * i->tag.doubleByOffset(((IntegratorScalar*) data)->m_force_offset[force_index]);

      );
}


void IntegratorScalar::integrateStep2()
{
}


#ifdef _OPENMP
string IntegratorScalar::dofIntegr() {
  return m_scalar_name;
}


void IntegratorScalar::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    p->tag.doubleByOffset(m_force_offset[force_index]) += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos];
    (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos] = 0;
  }
}

#endif
