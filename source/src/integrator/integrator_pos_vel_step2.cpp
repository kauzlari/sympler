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
#include "integrator_pos_vel_step2.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorPosVelStep2> integrator_pos_vel_step2("IntegratorPosVelStep2");

//---- Constructors/Destructor ----

IntegratorPosVelStep2::IntegratorPosVelStep2(Controller *controller):IntegratorPosition(controller)
{
  m_disp.assign(0);
  init();
}


IntegratorPosVelStep2::~IntegratorPosVelStep2()
{
}


//---- Methods ----

void IntegratorPosVelStep2::init()
{
  // some modules need to know whether there is an Integrator,
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorPosVelStep2");

  m_properties.setDescription
    ("Integrates the position r(t), velocity v(t), and displacement "
     "d(t) = r(t) - r(t = 0)  of each particle "
     "according to the following scheme:\n"
     "integration-step1: no activity\n"
     "integration-step2: v(t + dt) = v(t) + dt * F(t) / m\n"
     "                   r(t + dt) = r(t) + dt * v(t + dt)\n"
     "Here, F and m are particle force and mass, respectively, "
     "t is time and dt is the size of the integration "
     "time step (defined in the Controller). Further information on the "
     "integration-steps, including their place in the total SYMPLER "
     "workflow, can be found with the help option \"--help workflow\"."
     );

  STRINGPC
    (displacement, m_displacement_name,
     "Full name of the displacement, usable as attribute in other modules");
  
  STRINGPC
    (symbol, m_displacement_symbol,
     "Symbol assigned to the displacement, usable in algebraic expressions");
  
  m_displacement_name = "displacement";
  m_displacement_symbol = "ds";
  
}


void IntegratorPosVelStep2::setup()
{
  IntegratorPosition::setup();

  m_displacement_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      (m_displacement_name,
       DataFormat::POINT,
       true,
       m_displacement_symbol).offset;
}


void IntegratorPosVelStep2::isAboutToStart()
{
  IntegratorPosition::isAboutToStart();
}


void IntegratorPosVelStep2::integrateStep1()
{

}


void IntegratorPosVelStep2::integrateStep2()
{
  Phase *phase = M_PHASE;

  Controller* controller = M_CONTROLLER;
  
  size_t force_index = controller -> forceIndex();

  // velocity integration
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,

     i->v += i->dt * i->force[force_index] / i->m_mass;

     );

  // Currently(2018-06-27), this only invalidates the centre of mass
  // velocity
  phase -> invalidateVelocities();
  
  // position integration (will internally call integratePosition etc.)
  phase -> invalidatePositions(this);

  // Currently (2018-06-28) the Controller does not do this by itself
  // after integration step 2 (but after step 1), hence we do it here
  controller -> triggerNeighbourUpdate();
  
}


void IntegratorPosVelStep2::integratePosition(Particle* p, Cell* cell)
{
  size_t force_index;
  force_index = ((Controller*) m_parent) -> forceIndex();

  // is not 0. but should be irrelevant for the collision algorithm
  // since the new velocity is the only input for the new position
  point_t accel = { 0., 0., 0. };
  // point_t accel = p->force[force_index] / p->m_mass;

  // will also compute m_disp in IntegratorPosVelStep2::hitPos
  cell->doCollision(p, p->r, p->v, accel, (IntegratorPosition*) this);

  p->tag.pointByOffset(this->m_displacement_offset) += m_disp - p->r;
  
  // assuming the velocity is already the new one
  p->r += p->dt * p->v;

  p->tag.pointByOffset(this->m_displacement_offset) += p->r;
  // for next usage
  m_disp.assign(0);

  // inconsistency due to periodic BCs between the displacement and 
  // the position can not occur because the PBCs are only checked afterwards 
}


void IntegratorPosVelStep2::integrateVelocity(Particle* p)
{
  // Cell::updatePositions(IntegratorPosition* integrator) calls this
  // function AFTER IntegratorPosVelStep2::integratePosition(..), hence
  // too late for this algorithm. So we do the velocity update already
  // in IntegratorPosVelStep2::integrateStep2().
}


void IntegratorPosVelStep2::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results)
{
  double b, c;

  point_t surface_normal = wallTriangle->normal();

  b = surface_normal*p->v;
  c = surface_normal*p->r - wallTriangle->nDotR();

  // WallTriangle will take correct care of negative times
  results->push_back(-c/b);

}


void IntegratorPosVelStep2::hitPos
(const double& dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  m_disp = dt*p->v;

  hit_pos = p->r + m_disp;
}


#ifdef _OPENMP
string IntegratorPosVelStep2::dofIntegr() {
  return "vel_pos";
}

void IntegratorPosVelStep2::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->force[force_index][i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];

      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }

  }
}

#endif

