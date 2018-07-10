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
#include "controller.h"
#include "simulation.h"
#include "integrator_velocity_verlet_x_full.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()


const Integrator_Register<IntegratorVelocityVerletXFull> integrator_velocity_verlet_x_full("IntegratorVelocityVerletXFull");

//---- Constructors/Destructor ----

IntegratorVelocityVerletXFull::IntegratorVelocityVerletXFull(Controller *controller):IntegratorVelocityVerletX(controller)
{
  init();
}


IntegratorVelocityVerletXFull::~IntegratorVelocityVerletXFull()
{
}


//---- Methods ----

void IntegratorVelocityVerletXFull::init()
{
  // some modules need to know whether there is an Integrator, 
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorVelocityVerletXFull");

  m_properties.setDescription("Integrates the position and velocity coordinates of each particle according to the Velocity-Verlet Algorithm. The position is integrated by using a user-defined velocity. This is also the velocity, which is integrated.");

}

void IntegratorVelocityVerletXFull::setup()
{
  Integrator::setup();

  if(m_vSymbol == "undefined")
    throw gError("IntegratorVelocityVerletX::setup", "Attribute 'velocity' has value \"undefined\".");
  
  // the Integrator adds the symbol, because it is created before the other symbols;
  // the rest is the job of calculators and caches. If there are none, the velocity 
  // should be always zero
  if(!Particle::s_tag_format[m_colour].attrExists(m_vSymbol))
    /* !!! It IS an integrated quantity, so it IS persistent !!! 
    (notice the difference to IntegratorVelocityVerletX !!!)*/
    m_v_offset = 
        Particle::s_tag_format[m_colour].addAttribute(m_vSymbol, DataFormat::POINT,
      true).offset;
  else
    throw gError("IntegratorVelocityVerletX::setup", "Cannot add symbol " + m_vSymbol + " to species " + m_species + " because it is already used.");
  
}


void IntegratorVelocityVerletXFull::integrateStep2()
{
//   MSG_DEBUG("IntegratorVelocityVerletXFull::integrateStep2", "START: m_lambda = " << m_lambda << ", v_offset = " << m_v_offset);
  Phase *phase = M_PHASE;
  /*  m_force_index = M_CONTROLLER->forceIndex();*/
  size_t force_index = M_CONTROLLER->forceIndex();
  
//   m_other_force_index = (m_force_index+1)&(FORCE_HIST_SIZE-1);
  size_t other_force_index = (/*m_*/force_index+1)&(FORCE_HIST_SIZE-1);

  // if it was an estimate with m_lambda != 0.5 then make m_lambda = 0.5 now
  if (m_lambda != 0.5) {
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
        (phase, m_colour, this,
         i->tag.pointByOffset(m_v_offset) += 
             i->dt * /*((IntegratorPosition*) data)->*/m_lambda_diff * 
             i->force[/*((IntegratorPosition*) data)->m_*/other_force_index];
        );
  }

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
//        MSG_DEBUG("IntegratorVelocityVerletXFull::integrateStep2", name() << "MAIN: v BEFORE = " << i->tag.pointByOffset(m_v_offset));  
       i->tag.pointByOffset(m_v_offset)
           += i->dt/2 * i->force[/*((IntegratorPosition*) data)->m_*/force_index];
//        MSG_DEBUG("IntegratorPosition::integrateStep2", "MAIN: f = " << i->force[force_index]);
//        MSG_DEBUG("IntegratorVelocityVerlet::integrateStep2", name() << "MAIN: v AFTER = " << i->tag.pointByOffset(m_v_offset));  
      );

   phase->invalidateVelocities();
}

void IntegratorVelocityVerletXFull::integrateVelocity(Particle* p)
{
  int force_index;
  force_index = ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();
//   MSG_DEBUG("IntegratorVelocityVerlet::integrateVelocity", name() << "v BEFORE = " << p->v);  
  p->tag.pointByOffset(m_v_offset) += m_lambda * (p->dt * p->force[force_index]);
//   MSG_DEBUG("IntegratorVelocityVerlet::integrateVelocity", name() << "v AFTER = " << p->v);  
}
