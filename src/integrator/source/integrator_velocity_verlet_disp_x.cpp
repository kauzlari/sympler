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
#include "integrator_velocity_verlet_disp_x.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()


const Integrator_Register<IntegratorVelocityVerletDispX> integrator_velocity_verlet_disp_x("IntegratorVelocityVerletDispX");

//---- Constructors/Destructor ----

IntegratorVelocityVerletDispX::IntegratorVelocityVerletDispX(Controller *controller):IntegratorVelocityVerlet(controller)
{
  m_disp.assign(0);
  init();
}


IntegratorVelocityVerletDispX::~IntegratorVelocityVerletDispX()
{
}


//---- Methods ----

/*!
 * 
 */
void IntegratorVelocityVerletDispX::init()
{
  // some modules need to know whether there is an Integrator, 
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorVelocityVerletDispX");

  m_properties.setDescription("Integrates the position and momentum coordinates of each particle according to the Velocity-Verlet Algorithm and calculates the displacement for each particle. The velocity used for position and displacement integration is a user-defined symbol, instead of the usual integrated velocity.");


  STRINGPC
    (displacement, m_displacement_name,
     "Full name of the displacement, usable as attribute in other modules");
  
  STRINGPC
    (symbol, m_displacement_symbol,
     "Symbol assigned to the displacement, usable in algebraic expressions");

  m_displacement_name = "displacement";
  m_displacement_symbol = "ds";

  STRINGPC(velocity, m_vSymbol, "The symbol of the velocity used for integration of the positions.");
      
  m_vSymbol = "undefined";


}


void IntegratorVelocityVerletDispX::setup()
{
  Integrator::setup();

  if(m_vSymbol == "undefined")
    throw gError("IntegratorVelocityVerletDispX::setup", "Attribute 'velocity' has value \"undefined\".");
  
  // the Integrator adds the symbol, because it is created before the other symbols;
  // the rest is the job of calculators and caches. If there are none, the velocity 
  // should be always zero
  if(!Particle::s_tag_format[m_colour].attrExists(m_vSymbol))
    m_v_offset = 
      Particle::s_tag_format[m_colour].addAttribute(m_vSymbol, DataFormat::POINT,
      false /*it's not an integrated quantity, so not persistent !!!*/).offset;
  else
    throw gError("IntegratorVelocityVerletDispX::setup", "Cannot add symbol " + m_vSymbol + " to species " + m_species + " because it is already used.");
  
  m_displacement_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      (m_displacement_name,
       DataFormat::POINT,
       true,
       m_displacement_symbol).offset;
}


void IntegratorVelocityVerletDispX::integratePosition(Particle* p, Cell* cell)
{
  size_t force_index; 
  force_index = ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();
  
  point_t& vel = p->tag.pointByOffset(m_v_offset);
  point_t oldVel = vel;

  const point_t& a = p->force[force_index]/m_mass;

  cell->doCollision(p, p->r, vel, a, (IntegratorPosition*) this); 

  // Collision happened! So we have to change p->v! Currently we change it 
  // to vel, i.e., the user-defined velocity
  if(!(vel == oldVel))
    p->v = vel;

  p->tag.pointByOffset(this->m_displacement_offset) += m_disp - p->r;

  p->r += p->dt * (p->tag.pointByOffset(m_v_offset) + 0.5 * p->dt * a);  

  p->tag.pointByOffset(this->m_displacement_offset) += p->r;
  // for next usage
  m_disp.assign(0);

  // inconsistency due to periodic BCs between the displacement and 
  // the position can not occur because the PBCs are only checked afterwards

}

void IntegratorVelocityVerletDispX::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results)
{
  double a, b, c;
  double t0, t1;
  int n;
  point_t surface_normal = wallTriangle->normal();

  a = surface_normal*force / 2;
  b = surface_normal*p->tag.pointByOffset(m_v_offset);
  c = surface_normal*p->r - wallTriangle->nDotR();


  if (a!=0)
  { 
    n = gsl_poly_solve_quadratic(a, b, c, &t0, &t1);

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


void IntegratorVelocityVerletDispX::hitPos(/*WallTriangle* wallTriangle, */double dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  hit_pos = p->r + dt*(p->tag.pointByOffset(m_v_offset) + (dt/2)*(force/m_mass));

  m_disp = hit_pos - p->r;
}


#ifdef _OPENMP
string IntegratorVelocityVerletDispX::dofIntegr() {
  return "vel_pos";
}


#endif

