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


//
// C++ Implementation: integrator_rigidbodyMC
//
// Description: This integrator performs random moves of a rigid body for all its
// six degrees of freedom
//

#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_rigidbodyMC.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

//
#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorRigidBodyMC> integrator_rigidbodyMC("IntegratorRigidBodyMC");

//---- Constructors/Destructor ----

IntegratorRigidBodyMC::IntegratorRigidBodyMC(Controller *controller)
:Integrator(controller),m_first_sweep(true)
{
  init();
}


IntegratorRigidBodyMC::~IntegratorRigidBodyMC()
{
}


//---- Methods ----

void IntegratorRigidBodyMC::init()
{
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorRigidBodyMC");

  m_properties.setDescription("Integrates the the position and the orientation of selected rigid bodies according to an MC algorithm");
  DOUBLEPC
      (DeltaR, m_DeltaR,0,
       "Variation of the position (should be ca. 0.15 Angstrom).");
    m_DeltaR = 0.05;
  DOUBLEPC
      (PHIMax, m_PHIMax,0,
       "Maximum angle of rotation w.r.t. an arbitrary axis.");
    m_PHIMax = 0.1;
  DOUBLEPC
      (kT, m_kT,0,
       "Initial temperature.");
    m_kT = 1.0;
  STRINGPC
      (tensrotmat, m_rotmat_name,
       "Full name of the rotation matrix usable as attribute in other modules");
  STRINGPC
      (syrotmat, m_rotmat_symbol,
       "Symbol assigned to the rotation matrix, usable in algebraic expressions");
  m_rotmat_name = "rotmat";
  m_rotmat_symbol = "rotmat";

  STRINGPC
      (tensrotmatT, m_rotmatT_name,
       "Full name of the rotation matrix usable as attribute in other modules");
  STRINGPC
      (syrotmatT, m_rotmatT_symbol,
       "Symbol assigned to the rotation matrix, usable in algebraic expressions");
  m_rotmatT_name = "rotmatT";
  m_rotmatT_symbol = "rotmatT";

  STRINGPC
    (scq0, m_q0_name,
          "Full name of q0, usable as attribute in other modules.");
  STRINGPC
      (syq0, m_q0_symbol,
              "Symbol assigned to q0, usable in algebraic expressions.");
  m_q0_name = "q0";
  m_q0_symbol = "q0";

  STRINGPC
    (scq1, m_q1_name,
          "Full name of q1, usable as attribute in other modules.");
  STRINGPC
      (syq1, m_q1_symbol,
              "Symbol assigned to q1, usable in algebraic expressions.");
  m_q1_name = "q1";
  m_q1_symbol = "q1";

  STRINGPC
    (scq2, m_q2_name,
          "Full name of q2, usable as attribute in other modules.");
  STRINGPC
      (syq2, m_q2_symbol,
              "Symbol assigned to q2, usable in algebraic expressions.");
  m_q2_name = "q2";
  m_q2_symbol = "q2";

  STRINGPC
    (scq3, m_q3_name,
          "Full name of q3, usable as attribute in other modules.");
  STRINGPC
      (syq3, m_q3_symbol,
              "Symbol assigned to q3, usable in algebraic expressions.");
  m_q3_name = "q3";
  m_q3_symbol = "q3";
  STRINGPC
    (Etot,m_Etot_name,
     "Full name of total energy, usable as attribute in other modules.");
  m_Etot_name="undefined";
  STRINGPC
    (displacement, m_displacement_name,
     "Full name of the displacement, usable as attribute in other modules");
  
  STRINGPC
    (symbol, m_displacement_symbol,
     "Symbol assigned to the displacement, usable in algebraic expressions");

  m_displacement_name = "displacement";
  m_displacement_symbol = "ds";

}

void IntegratorRigidBodyMC::setup()
{
  M_CONTROLLER->registerForSetupAfterParticleCreation(this);
  Integrator::setup();

  DataFormat::attribute_t tmpAttr;

/*!
 * Setup the quaternion components
 */
  m_q0_offset =
    Particle::s_tag_format[m_colour].addAttribute
      (m_q0_name,
       DataFormat::DOUBLE,
       true,
       m_q0_symbol).offset;
  m_q1_offset =
    Particle::s_tag_format[m_colour].addAttribute
      (m_q1_name,
       DataFormat::DOUBLE,
       true,
       m_q1_symbol).offset;
  m_q2_offset =
    Particle::s_tag_format[m_colour].addAttribute
      (m_q2_name,
       DataFormat::DOUBLE,
       true,
       m_q2_symbol).offset;
  m_q3_offset =
    Particle::s_tag_format[m_colour].addAttribute
      (m_q3_name,
       DataFormat::DOUBLE,
       true,
       m_q3_symbol).offset;
/*!
 * Setup the rotation matrix
 */
  m_rotmat_offset =
      Particle::s_tag_format[m_colour].addAttribute
      (m_rotmat_name,
       DataFormat::TENSOR,
       true,
       m_rotmat_symbol).offset;
/*!
 * Setup the transpose of rotation matrix
 */
  m_rotmatT_offset =
      Particle::s_tag_format[m_colour].addAttribute
      (m_rotmatT_name,
       DataFormat::TENSOR,
       true,
       m_rotmatT_symbol).offset;
/*!
 * Setup the old position
 */
  m_oldpos_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      (m_oldpos_name,
       DataFormat::POINT,
       true,m_oldpos_symbol).offset;
/*!
 * Setup the displacement
 */
  m_displacement_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      (m_displacement_name,
       DataFormat::POINT,
       true,
       m_displacement_symbol).offset;
/*!
 * Fetch the timestep from controller
 */
  m_dt = M_CONTROLLER->dt();
}
/*!
 *  Calculate rotation matrix after particle setup
 */
void IntegratorRigidBodyMC::setupAfterParticleCreation()
{
  m_reject = false;
  EtotOld = 1.0e40;
/*!
 * Fetch the total energy from the particleScalar
 * Why is this here? ::setup is too early we must
 * wait because integrators are setup before all 
 * other modules.
 */
  if(m_Etot_name == "undefined")
    throw gError("IntegratorRigidBodyMC::setup", "Attribute 'Etot' has value \"undefined\". Please give the name under which the displacement will be computed by another module, e.g., by a ParticleScalar.");

  if(Particle::s_tag_format[m_colour].attrExists(m_Etot_name)) {
	DataFormat::attribute_t attrEtot = Particle::s_tag_format[m_colour].attrByName(m_Etot_name);
	
	if(attrEtot.datatype != DataFormat::DOUBLE)
	  throw gError("IntegratorRigidBodyMC::setup", "the symbol " + m_Etot_name + " is not registered as a scalar.");
	
	m_Etot_o = attrEtot.offset;
      }
      else {
	throw gError("IntegratorRigidBodyMC::setup", "No symbol '" + m_Etot_name + "' found for species '" + m_species + "'! You need another module introducing it. This might be for example a ParticleScalar.");
      }
// Let us define some helpers for writing the formulae
// 1. the 4 components of the quaternion 
//
#define Q0 (i->tag.doubleByOffset(m_q0_offset))
#define Q1 (i->tag.doubleByOffset(m_q1_offset))
#define Q2 (i->tag.doubleByOffset(m_q2_offset))
#define Q3 (i->tag.doubleByOffset(m_q3_offset))
#define  X (i->r.x)
#define  Y (i->r.y)
#define  Z (i->r.z)
//
//
// 7. We define the rotation matrix here
//
#define R (i->tag.tensorByOffset(m_rotmat_offset))
#define RT (i->tag.tensorByOffset(m_rotmatT_offset))
//
  Phase *phase = M_PHASE;
  FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, m_colour, this,
      R(0,0) = Q0*Q0+Q1*Q1-Q2*Q2-Q3*Q3;
      R(0,1) = 2*(Q1*Q2-Q0*Q3);
      R(0,2) = 2*(Q1*Q3+Q0*Q2);
      R(1,0) = 2*(Q1*Q2+Q0*Q3);
      R(1,1) = Q0*Q0-Q1*Q1+Q2*Q2-Q3*Q3;
      R(1,2) = 2*(Q2*Q3-Q0*Q1);
      R(2,0) = 2*(Q1*Q3-Q0*Q2);
      R(2,1) = 2*(Q2*Q3+Q0*Q1);
      R(2,2) = Q0*Q0-Q1*Q1-Q2*Q2+Q3*Q3;
// Here we define the transpose of R(i,j)
      RT(0,0) = R(0,0);
      RT(1,0) = R(0,1);
      RT(2,0) = R(0,2);
      RT(0,1) = R(1,0);
      RT(1,1) = R(1,1);
      RT(2,1) = R(1,2);
      RT(0,2) = R(2,0);
      RT(1,2) = R(2,1);
      RT(2,2) = R(2,2);

      );
  } 

/*!
 *
 */
void IntegratorRigidBodyMC::integrateStep1(){
  Phase *phase = M_PHASE;
  double PHI, DISP, norm, u, v, w, DQ0, DQ1, DQ2, DQ3, P0, P1, P2, P3;
  double rnumber, expfact;
//  
  EtotNew = 0.0;
  FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, m_colour, this,
// Here we calculate first the energy  
      EtotNew+=i->tag.doubleByOffset(((IntegratorRigidBodyMC*) data)->m_Etot_o);
      );
  EtotNew = EtotNew/M_PHASE->returnNofPartC(m_colour);
// Based on this new energy we decide if to reject or to accept
  if(EtotNew<EtotOld){
    m_reject = false;
    EtotOld = EtotNew;
  }
  else{
    rnumber=m_rng.uniform();
    expfact=exp(-(EtotNew-EtotOld)/m_kT);
    if( rnumber - expfact < 0.0){
      m_reject = false;
//MSG_DEBUG("IntegratorRigidBodyMC::integrateStep1()", name() << " r and exp " << rnumber << ", " << expfact);
    }
    else{
      m_reject = true;
//MSG_DEBUG("IntegratorRigidBodyMC::integrateStep1()", name() <<  " r and exp " << rnumber << ", " << expfact);
    }
  }

// 
  m_velmax = m_DeltaR/m_dt;
  if(m_reject){
    M_PHASE->invalidatePositions((IntegratorPosition*) this);
    m_reject = false;
  }
  M_PHASE->invalidatePositions((IntegratorPosition*) this);
//
//
  FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, m_colour, this,
// Here we first draw a random axis o rotate about
      u = 1.0-2.0*m_rng.uniform();
      v = 1.0-2.0*m_rng.uniform();
      w = 1.0-2.0*m_rng.uniform();
      norm = sqrt(u*u+v*v+w*w);
      u = u/norm;
      v = v/norm;
      w = w/norm;
// Then we draw the angle to turn the rigid body
      PHI = m_PHIMax*(.5-m_rng.uniform());
// Carful this is already phi/2 accroding to definition of 
// quaternions: q0 = cos(phi/2), q1 = sin(phi/2)*u, ...
      DQ0 = cos(PHI);
      DQ1 = sin(PHI)*u;
      DQ2 = sin(PHI)*v;
      DQ3 = sin(PHI)*w;
      norm = sqrt(DQ0*DQ0+DQ1*DQ1+DQ2*DQ2+DQ3*DQ3);
      DQ0 = DQ0/norm;
      DQ1 = DQ1/norm;
      DQ2 = DQ2/norm;
      DQ3 = DQ3/norm;
// And now we calculate the new quaternion Q' by rotation
// with the DQi. Formula Q'=DQ*Q*DQ^-1 (* = quaternion multiplication)
// First the DQ*Q = P
      P0 = DQ0*Q0-DQ1*Q1-DQ2*Q2-DQ3*Q3;
      P1 = DQ0*Q1+DQ1*Q0+DQ2*Q3-DQ3*Q2;
      P2 = DQ0*Q2+DQ2*Q0+DQ3*Q1-DQ1*Q3;
      P3 = DQ0*Q3+DQ3*Q0+DQ1*Q2-DQ2*Q1;
// Then the P*DQ^-1 
      Q0 =  P0*Q0+P1*Q1+P2*Q2+P3*Q3;
      Q1 = -P0*Q1+P1*Q0-(P2*Q3-P3*Q2);
      Q2 = -P0*Q2+P2*Q0-(P3*Q1-P1*Q3);
      Q3 = -P0*Q3+P3*Q0-(P1*Q2-P2*Q1);
      norm = sqrt(Q0*Q0+Q1*Q1+Q2*Q2+Q3*Q3);
      Q0 = Q0/norm;
      Q1 = Q1/norm;
      Q2 = Q2/norm;
      Q3 = Q3/norm;
//  
      R(0,0) = Q0*Q0+Q1*Q1-Q2*Q2-Q3*Q3;
      R(0,1) = 2*(Q1*Q2-Q0*Q3);
      R(0,2) = 2*(Q1*Q3+Q0*Q2);
      R(1,0) = 2*(Q1*Q2+Q0*Q3);
      R(1,1) = Q0*Q0-Q1*Q1+Q2*Q2-Q3*Q3;
      R(1,2) = 2*(Q2*Q3-Q0*Q1);
      R(2,0) = 2*(Q1*Q3-Q0*Q2);
      R(2,1) = 2*(Q2*Q3+Q0*Q1);
      R(2,2) = Q0*Q0-Q1*Q1-Q2*Q2+Q3*Q3;
      // Here we define the transpose of R(i,j)
      RT(0,0) = R(0,0);
      RT(1,0) = R(0,1);
      RT(2,0) = R(0,2);
      RT(0,1) = R(1,0);
      RT(1,1) = R(1,1);
      RT(2,1) = R(1,2);
      RT(0,2) = R(2,0);
      RT(1,2) = R(2,1);
      RT(2,2) = R(2,2);
      );
}

void IntegratorRigidBodyMC::integratePosition(Particle* p, Cell* cell)
{
  double randvel, u, v, w, norm;
// Currently (2010-05-05), pt is a const point& argument, so using it in the p->r += ... line is safe
  const point_t pt = {0,0,0}; 
// Here we displace the body by a random vector 
// we first draw a random unit vector as above
  if(m_reject){
    p->tag.pointByOffset(m_displacement_offset)=p->tag.pointByOffset(m_oldpos_offset)-p->r;
    p->r=p->tag.pointByOffset(m_oldpos_offset);
  }
  else{
      p->tag.pointByOffset(m_oldpos_offset)=p->r;
      u = 1.0-2.0*m_rng.uniform();
      v = 1.0-2.0*m_rng.uniform();
      w = 1.0-2.0*m_rng.uniform();
      norm = sqrt(u*u+v*v+w*w);
      u = u/norm;
      v = v/norm;
      w = w/norm;
      randvel = m_velmax*m_rng.uniform();
      p->v.x = randvel*u;
      p->v.y = randvel*v;
      p->v.z = randvel*w;
  cell->doCollision(p, p->r, p->v, pt, (IntegratorPosition*) this);
  p->r += p->dt * p->v ;
  }
// MSG_DEBUG("IntegratorRigidBodyMC::integratePosition", name() << "doCollision called, final update of r done");
}

void IntegratorRigidBodyMC::integrateStep2(){
}

void IntegratorRigidBodyMC::solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t &force, vector<double>* results)
{
  double  b, c;
  double t0, t1;
  point_t surface_normal = wallTriangle->normal();

  b = surface_normal*p->v;
  c = surface_normal*p->r - wallTriangle->nDotR();
  t0 = -c/b;
  results->push_back(t0);
}


void IntegratorRigidBodyMC::hitPos(/*WallTriangle* wallTriangle, */double dt, const Particle* p, point_t &hit_pos, const point_t &force)
{
  hit_pos = p->r + dt*p->v;
}

#ifdef _OPENMP
string IntegratorRigidBodyMC::dofIntegr() {
  return m_omega_name;
}

void IntegratorRigidBodyMC::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->tag.pointByOffset(m_omega_force_offset[force_index])[i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
  }
}

#endif

