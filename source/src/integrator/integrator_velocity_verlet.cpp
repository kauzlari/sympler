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
#include "integrator_velocity_verlet.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorVelocityVerlet> integrator_velocity_verlet("IntegratorVelocityVerlet");

//---- Constructors/Destructor ----

IntegratorVelocityVerlet::IntegratorVelocityVerlet(Controller *controller):IntegratorPosition(controller)
{
  init();
}


IntegratorVelocityVerlet::~IntegratorVelocityVerlet()
{
}


//---- Methods ----

void IntegratorVelocityVerlet::init()
{
  // some modules need to know whether there is an Integrator,
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorVelocityVerlet");

  m_properties.setDescription("Integrates the position and momentum coordinates of each particle according to the Velocity-Verlet Algorithm");

    DOUBLEPC
    (lambda,
     m_lambda,
     0,
     "Lambda parameter for modified velocity verlet algorithm.");

  m_lambda = 0.5;

}


void IntegratorVelocityVerlet::isAboutToStart()
{
  IntegratorPosition::isAboutToStart();

  m_lambda_diff = 0.5 - m_lambda;
}


void IntegratorVelocityVerlet::integrateStep1()
{
  M_PHASE->invalidatePositions((IntegratorPosition*) this);
}


void IntegratorVelocityVerlet::integrateStep2()
{
  Phase *phase = M_PHASE;
/*  m_force_index = M_CONTROLLER->forceIndex();*/
  size_t force_index = M_CONTROLLER->forceIndex();

//   m_other_force_index = (m_force_index+1)&(FORCE_HIST_SIZE-1);
  size_t other_force_index = (/*m_*/force_index+1)&(FORCE_HIST_SIZE-1);

  // if it was an estimate with m_lambda != 0.5 then make m_lambda = 0.5 now
  if (m_lambda != 0.5) {
    FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->v += i->dt * /*((IntegratorPosition*) data)->*/m_lambda_diff *
           i->force[/*((IntegratorPosition*) data)->m_*/other_force_index]/m_mass;
      );
  }

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
// if (i->mySlot==35)
//      MSG_DEBUG("IntegratorVelocityVerlet::integrateStep2", name() << "v BEFORE = " << i->force[force_index]);
     i->v += i->dt/2 * i->force[/*((IntegratorPosition*) data)->m_*/force_index]/m_mass;
// if (i->mySlot == 0){
// MSG_DEBUG("IntegratorVelocityVerlet::integrateStep2", "velocity = " << i->v);
// MSG_DEBUG("IntegratorVelocityVerlet::integrateStep2", "force = " << i->force[/*((IntegratorPosition*) data)->m_*/force_index]);}
//      MSG_DEBUG("IntegratorPosition::integrateStep2", "f = " /*<< i->force[((IntegratorPosition*) data)->m_force_index]*/);
// if (i->mySlot==35)
//      MSG_DEBUG("IntegratorVelocityVerlet::integrateStep2", name() << "v AFTER = " << i->force[force_index]);
    );

   phase->invalidateVelocities();
}


void IntegratorVelocityVerlet::integratePosition(Particle* p, Cell* cell)
{
  size_t force_index;
  force_index = ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();

  point_t accel = p->force[force_index]/m_mass;
  // Currently (2010-05-05), pt is a const point& argument, so using it in the p->r += ... line is safe
  cell->doCollision(p, p->r, p->v, accel, (IntegratorPosition*) this);

  p->r += p->dt * (p->v + 0.5 * p->dt * accel);
  //MSG_DEBUG("IntegratorVelocityVerlet::integratePosition", name() << "pos_force= " <<  p->force[force_index]);
  //MSG_DEBUG("IntegratorVelocityVerlet::integratePosition", name() << "position= " <<  p->r);

}


void IntegratorVelocityVerlet::integrateVelocity(Particle* p)
{
  size_t force_index;
  // FIXME: inefficient to set the force_index again and again for each
  // particle. Can we use Node::precompute() or s.th. similar?
  force_index = ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();
//   MSG_DEBUG("IntegratorVelocityVerlet::integrateVelocity", name() << "v BEFORE = " << p->v);
  p->v += m_lambda * (p->dt * p->force[force_index]/m_mass);
  //MSG_DEBUG("IntegratorVelocityVerlet::integrateVelocity", name() << "v AFTER = " << p->v);
}


void IntegratorVelocityVerlet::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results)
{
  double a, b, c;
  double t0, t1;
  int n;
  point_t surface_normal = wallTriangle->normal();

  a = surface_normal*force/m_mass / 2;
  b = surface_normal*p->v;
  c = surface_normal*p->r - wallTriangle->nDotR();


  if (a!=0)
  {
    n = gsl_poly_solve_quadratic(a, b, c, &t0, &t1);

//     MSG_DEBUG("IntegratorVelocityVerlet::solveHitTimeEquation", "n = " << n << ", t0 = " << t0 << ", t1 = " << t1 << " for a = " << a << ", b = " << b << ", c = " << c << ", force = " << force << ", v = " << p->v << ", r = " << p->r << ", ndotr = " << wallTriangle->nDotR() << ", surfnormal = " << surface_normal);

    if (n == 0 || (t0 < c_wt_time_eps && t1 < c_wt_time_eps))
    {
    }

    else
    {
      if (t0 < c_wt_time_eps)
      {
        t0 = t1;
        results->push_back(t0);
      }

    results->push_back(t0);
    results->push_back(t1);
    sort(results->begin(), results->end());
    }
  }
  else
  {
    t0 = -c/b;

    if (t0 < c_wt_time_eps)
    {
      n = 0;
    }

    n = 1;
    results->push_back(t0);
  }
}


void IntegratorVelocityVerlet::hitPos(/*WallTriangle* wallTriangle, */double dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  hit_pos = p->r + dt*p->v + dt*dt/2*force/m_mass;
}


#ifdef _OPENMP
string IntegratorVelocityVerlet::dofIntegr() {
  return "vel_pos";
}

void IntegratorVelocityVerlet::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->force[force_index][i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
// MSG_DEBUG("IntegratorVelocityVerlet::mergeCopies", " real force to be added = " << (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] << " slot = " << p->mySlot);
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
//    MSG_DEBUG("IntegratorVelocityVerlet::mergeCopies", " force after merge = " << p->force[force_index]);
  }
}

#endif

