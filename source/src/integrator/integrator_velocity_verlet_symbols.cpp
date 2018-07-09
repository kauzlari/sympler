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




#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_velocity_verlet_symbols.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorVelocityVerletSymbols> integrator_velocity_verlet_symbols("IntegratorVelocityVerletSymbols");

//---- Constructors/Destructor ----

IntegratorVelocityVerletSymbols::IntegratorVelocityVerletSymbols(Controller *controller):Integrator(controller)
{
  init();
}


IntegratorVelocityVerletSymbols::~IntegratorVelocityVerletSymbols()
{
}


//---- Methods ----

void IntegratorVelocityVerletSymbols::init()
{
  m_properties.setClassName("IntegratorVelocityVerletSymbols");
  m_properties.setName("IntegratorVelocityVerletSymbols");

  m_properties.setDescription("Integrates \"positions\" and \"momenta\" stored as symbols in each particle's tag according to the Velocity-Verlet Algorithm. Hence, the true positions and velocities are not affected and the particle does NOT move. ");

  DOUBLEPC
    (lambda,
     m_lambda,
     0,
     "Lambda parameter for modified velocity verlet algorithm.");

  m_lambda = 0.5;

  DOUBLEPC
    (mass, m_mass,0,
     "Particle mass of the species this integrator is intended for. Default mass = 1. Pay attention the mass is only "
     "effective for the integrator and does not affect the thermostat for instance");
  m_mass = 1;

  STRINGPC
      (velName, m_vel_name,
       "Full name of the velocity, usable as attribute in other modules");

  STRINGPC
      (velSymbol, m_vel_symbol,
       "Symbol name for the velocity, usable in algebraic expressions");

  m_vel_name = "velocity";
  m_vel_symbol = "vel";

  STRINGPC
      (posName, m_pos_name,
       "Full name of the position, usable as attribute in other modules");

  STRINGPC
      (posSymbol, m_pos_symbol,
       "Symbol name for the position, usable in algebraic expressions");

  m_pos_name = "velocity";
  m_pos_symbol = "vel";

}


void IntegratorVelocityVerletSymbols::setup()
{
  Integrator::setup();

  DataFormat::attribute_t tmpAttr;

  m_pos_offset =
      Particle::s_tag_format[m_colour].addAttribute
      (m_pos_name,
       DataFormat::POINT,
       true,
       m_pos_symbol).offset;

  m_vel_offset =
      Particle::s_tag_format[m_colour].addAttribute
      (m_vel_name,
       DataFormat::POINT,
       true,
       m_vel_symbol).offset;

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    tmpAttr =
        Particle::s_tag_format[m_colour].addAttribute
        (STR_FORCE + STR_DELIMITER + m_vel_name + STR_DELIMITER + ObjToString(i),
         DataFormat::POINT,
         true);

    m_force_offset[i] = tmpAttr.offset;
    m_fAttr_index[i] = tmpAttr.index;
  }
  m_dt = M_CONTROLLER->dt();

}

void IntegratorVelocityVerletSymbols::isAboutToStart()
{

  Phase *phase = M_PHASE;
  double invMass;

  if(m_mass <=0)
    throw gError("IntegratorVelocityVerletSymbols::isAboutToStart", "Invalid value \"" + ObjToString(m_mass) + "\" for attribute 'mass'. Must be >0!");
  invMass = 1/m_mass;
  m_dt_div_mass = m_dt*invMass;
  m_dt_div2_mass = m_dt_div_mass / 2;
  m_lambda_diff = 0.5 - m_lambda;

  MSG_DEBUG("IntegratorVelocityVerletSymbols::isAboutToStart", "m_colour = " << m_colour);

  size_t counter = 0;
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       ++counter;
       for (int j = 0; j < FORCE_HIST_SIZE; ++j)
           i->tag.pointByOffset(m_force_offset[j]).assign(0);
      );
  if(counter == 0)
    throw gError("IntegratorVelocityVerletSymbols::isAboutToStart", "no free particles found for species " + m_species + "! Don't instantiate an Integrator for positions and velocities in that case. Use another module to create the species.");
  // FIXME: so we need some SpeciesCreator to make it more transparent

  // FIXME: put all in this function into the general setup for Nodes after the particle creation or into s.th. even more general

}

void IntegratorVelocityVerletSymbols::unprotect(size_t index)
{
  Phase *phase = M_PHASE;

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->tag.unprotect(m_fAttr_index[index]);
       if(!index) i->tag.protect(m_fAttr_index[FORCE_HIST_SIZE-1]);
       else i->tag.protect(m_fAttr_index[index-1]);
      );
}

void IntegratorVelocityVerletSymbols::integrateStep1()
{
  Phase *phase = M_PHASE;
  size_t force_index = ((Controller*) m_parent)->forceIndex();

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     // positions
     i->tag.pointByOffset(m_pos_offset) += m_dt * (i->tag.pointByOffset(m_vel_offset) + m_dt_div2_mass * i->tag.pointByOffset(m_force_offset[force_index]));
  
     // velocities: prediction 
     i->tag.pointByOffset(m_vel_offset) += m_lambda * (m_dt_div_mass * i->tag.pointByOffset(m_force_offset[force_index]));

     );
}


void IntegratorVelocityVerletSymbols::integrateStep2()
{
  Phase *phase = M_PHASE;
  size_t force_index = M_CONTROLLER->forceIndex();
  size_t other_force_index = (force_index+1)&(FORCE_HIST_SIZE-1);

  // if it was an estimate with m_lambda != 0.5 then make m_lambda = 0.5 now
  if (m_lambda != 0.5) {
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->tag.pointByOffset(m_vel_offset) += m_lambda_diff * m_dt_div_mass
       * i->tag.pointByOffset(m_force_offset[other_force_index]);
      );
  }

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     i->tag.pointByOffset(m_vel_offset) += m_dt_div2_mass * i->tag.pointByOffset(m_force_offset[force_index]);
    );

}

#ifdef _OPENMP
string IntegratorVelocityVerletSymbols::dofIntegr() {
  // returns the name of the velocity
  return m_vel_name;
}

void IntegratorVelocityVerletSymbols::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->tag.pointByOffset(m_force_offset[force_index])[i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
// MSG_DEBUG("IntegratorVelocityVerletSymbols::mergeCopies", " real force to be added = " << (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] << " slot = " << p->mySlot);
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
//    MSG_DEBUG("IntegratorVelocityVerletSymbols::mergeCopies", " force after merge = " << p->force[force_index]);
  }
}


#endif

