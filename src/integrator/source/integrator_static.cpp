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
#include "integrator_static.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()


const Integrator_Register<IntegratorStatic> integrator_static("IntegratorStatic");

//---- Constructors/Destructor ----

IntegratorStatic::IntegratorStatic(Controller *controller): IntegratorPosition(controller)
{
    init();
}


IntegratorStatic::~IntegratorStatic()
{
}



//---- Methods ----

void IntegratorStatic::init()
{
  m_properties.setClassName("IntegratorPosition");

  m_properties.setName("IntegratorStatic");

  m_properties.setDescription(
    "Integrates the position coordinates of each particle according to the user-defined particle velocity.."
  );


  FUNCTIONFIXEDPC(velX, m_velocity_x,
                  "The x-component of the velocity. You can include time dependency with the variable 't' and spatial dependency with 'x', 'y', 'z'.");

  m_velocity_x.addVariable("x");
  m_velocity_x.addVariable("y");
  m_velocity_x.addVariable("z");
  m_velocity_x.addVariable("t");
  m_velocity_x.setExpression("0");

  FUNCTIONFIXEDPC(velY, m_velocity_y,
                  "The x-component of the velocity. You can include time dependency with the variable 't' and spatial dependency with 'x', 'y', 'z'.");

  m_velocity_y.addVariable("x");
  m_velocity_y.addVariable("y");
  m_velocity_y.addVariable("z");
  m_velocity_y.addVariable("t");
  m_velocity_y.setExpression("0");

  FUNCTIONFIXEDPC(velZ, m_velocity_z,
                  "The x-component of the velocity. You can include time dependency with the variable 't' and spatial dependency with 'x', 'y', 'z'.");

  m_velocity_z.addVariable("x");
  m_velocity_z.addVariable("y");
  m_velocity_z.addVariable("z");
  m_velocity_z.addVariable("t");
  m_velocity_z.setExpression("0");
}


void IntegratorStatic::isAboutToStart()
{
  m_dt = M_CONTROLLER->dt();
}


void IntegratorStatic::integrateStep1()
{
  M_PHASE->invalidatePositions((IntegratorPosition*) this);
}


void IntegratorStatic::integrateStep2()
{
}


void IntegratorStatic::integratePosition(Particle* p, Cell* cell)
{
  double time;
  time = ((Controller*) m_parent)->time();

  point_t& pos = p->r;
  point_t& vel = p->v;

  vel.x = m_velocity_x(pos.x, pos.y, pos.z, time);
  vel.y = m_velocity_y(pos.x, pos.y, pos.z, time);
  vel.z = m_velocity_z(pos.x, pos.y, pos.z, time);

  pos += p->dt * vel;
}


void IntegratorStatic::integrateVelocity(Particle* p)
{
}


void IntegratorStatic::setup()
{
  Integrator::setup();

}


void IntegratorStatic::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t
&force, vector<double>* results)
{
//   results->clear();
}


void IntegratorStatic::hitPos
(const double& dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
}


#ifdef _OPENMP
string IntegratorStatic::dofIntegr() {
  return "vel_pos";
}


void IntegratorStatic::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->force[force_index][i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
  }
}

#endif






