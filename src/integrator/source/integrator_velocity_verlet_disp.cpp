/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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
#include "integrator_velocity_verlet_disp.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()


const Integrator_Register<IntegratorVelocityVerletDisp> integrator_velocity_verlet_disp("IntegratorVelocityVerletDisp");

//---- Constructors/Destructor ----

IntegratorVelocityVerletDisp::IntegratorVelocityVerletDisp(Controller *controller):IntegratorVelocityVerlet(controller)
{
  m_disp.assign(0);
  init();
}


IntegratorVelocityVerletDisp::~IntegratorVelocityVerletDisp()
{
}


//---- Methods ----

/*!
 * 
 */
void IntegratorVelocityVerletDisp::init()
{
  // some modules need to know whether there is an Integrator, 
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorVelocityVerletDisp");

  m_properties.setDescription("Integrates the position and momentum coordinates of each particle according to the Velocity-Verlet Algorithm and calculates the displacement for each particle.");


  STRINGPC
    (displacement, m_displacement_name,
     "Full name of the displacement, usable as attribute in other modules");
  
  STRINGPC
    (symbol, m_displacement_symbol,
     "Symbol assigned to the displacement, usable in algebraic expressions");

  m_displacement_name = "displacement";
  m_displacement_symbol = "ds";

}


void IntegratorVelocityVerletDisp::setup()
{
  IntegratorVelocityVerlet::setup();

  m_displacement_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      (m_displacement_name,
       DataFormat::POINT,
       true,
       m_displacement_symbol).offset;
}


void IntegratorVelocityVerletDisp::integratePosition(Particle* p, Cell* cell)
{
  size_t force_index; 
  force_index = ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();

  point_t accel = p->force[force_index]/m_mass;

  // will also compute m_disp in IntegratorVelocityVerletDisp::hitPos
  cell->doCollision(p, p->r, p->v, accel, (IntegratorPosition*) this); 

  p->tag.pointByOffset(this->m_displacement_offset) += m_disp - p->r;

  p->r += p->dt * (p->v + 0.5 * p->dt * accel);  

  p->tag.pointByOffset(this->m_displacement_offset) += p->r;
  // for next usage
  m_disp.assign(0);

  // inconsistency due to periodic BCs between the displacement and 
  // the position can not occur because the PBCs are only checked afterwards

}


void IntegratorVelocityVerletDisp::hitPos
(const double& dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  m_disp = dt*(p->v + dt/2*force/m_mass);

  hit_pos = p->r + m_disp;
}


#ifdef _OPENMP
string IntegratorVelocityVerletDisp::dofIntegr() {
  return "vel_pos";
}


#endif

