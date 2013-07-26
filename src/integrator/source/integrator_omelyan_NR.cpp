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

//
// C++ Implementation: integrator_omelyan_NR
//
// Description: This integrator is based on the paper by Igor P. Omelyan,
// On the numerical integration of rigid polyatomics: The modified quaternion approach,
// Computers in Physics, Vol. 12, No. 1, p. 97 (1998).
// Careful, I use the notation Q=(q_0,\vect{q}) here as in Allen/Tidesley, Computer Simulation 
// of Liquids (Clarendon, Oxford, 1987). The difference to its simple version is a NewtonRaphson iteration
//

#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_omelyan_NR.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

//
#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorOmelyanNR> integrator_omelyan_NR("IntegratorOmelyanNR");

//---- Constructors/Destructor ----

IntegratorOmelyanNR::IntegratorOmelyanNR(Controller *controller)
:Integrator(controller),m_first_sweep(true)
{
  init();
}


IntegratorOmelyanNR::~IntegratorOmelyanNR()
{
}


//---- Methods ----

void IntegratorOmelyanNR::init()
{
  // some modules need to know whether there is an Integrator,
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorOmelyanNR");
  m_properties.setName("IntegratorOmelyanNR");

  m_properties.setDescription("Integrates the orientation and angular velocity of each rigid body attached to a particle according to the algorithm described by Omelyan");
  DOUBLEPC
  	    (J1, m_J1,0,
  	     "1st principle component of the inertial tensor of the species this integrator is intended for.");
  m_J1 = 1;
  DOUBLEPC
  	    (J2, m_J2,0,
  	     "2nd principle component of the inertial tensor of the species this integrator is intended for.");
  m_J2 = 1;
  DOUBLEPC
  	    (J3, m_J3,0,
  	     "3rd principle component of the inertial tensor of the species this integrator is intended for.");
  m_J3 = 1;

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
    (scomega, m_omega_name,
     "Full name of omega, usable as attribute in other modules.");
  STRINGPC
    (syomega, m_omega_symbol,
     "Symbol assigned to omega, usable in algebraic expressions.");
  m_omega_name = "omega";
  m_omega_symbol = "omega";
  BOOLPC
        (normalizeQuaternion, m_normalize,
	      "Normalize quaternions at the end of Step2.");
  m_normalize = false;
}

void IntegratorOmelyanNR::setup()
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
 * Setup the angular velocity components
 */
  m_omega_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      (m_omega_name,
       DataFormat::POINT,
       true,
       m_omega_symbol).offset;
/*!
 * Setup the torque components
 */
  m_torque_offset = 
    Particle::s_tag_format[m_colour].addAttribute
      ("m_torque",
       DataFormat::POINT,
       true,
       "m_torque").offset;
/*!
 * Setup the forces
 */
  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    /*!
     * The angular velocity forces setup
     */
    tmpAttr =
      Particle::s_tag_format[m_colour].addAttribute
       (STR_FORCE + STR_DELIMITER + m_omega_name + STR_DELIMITER + ObjToString(i),
       DataFormat::POINT,
       true);
    m_omega_force_offset[i] = tmpAttr.offset;
    m_omega_fAttr_index[i] = tmpAttr.index;
  }
/*!
 * Fetch the timestep from controller
 */
  m_dt = M_CONTROLLER->dt();
}
/*!
 *  Calculate rotation matrix after particle setup
 */
void IntegratorOmelyanNR::setupAfterParticleCreation()
{
// Let us define some helpers for writing the formulae
// 1. the 4 components of the quaternion 
//
#define Q0 (i->tag.doubleByOffset(m_q0_offset))
#define Q1 (i->tag.doubleByOffset(m_q1_offset))
#define Q2 (i->tag.doubleByOffset(m_q2_offset))
#define Q3 (i->tag.doubleByOffset(m_q3_offset))
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
 * Unprotect the forces on omega
 */
void IntegratorOmelyanNR::unprotect(size_t index)
{
  Phase *phase = M_PHASE;

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->tag.unprotect(m_omega_fAttr_index[index]);
       if(!index) i->tag.protect(m_omega_fAttr_index[FORCE_HIST_SIZE-1]);
       else i->tag.protect(m_omega_fAttr_index[index-1]);

      );
}

/*!
 *
 */
void IntegratorOmelyanNR::integrateStep1(){
  Phase *phase = M_PHASE;
  size_t force_index = M_CONTROLLER->forceIndex();
  size_t other_force_index = (force_index+1)&(FORCE_HIST_SIZE-1);
  double LAMBDA;
//
// Let us define some helpers for writing the formulae
//
// 2. the 3 components of the angular velocity
//
#define Omega (i->tag.pointByOffset(m_omega_offset))
//
//
// 3. the 3 components of the torque at actual time
//
#define TauOld (i->tag.pointByOffset(m_omega_force_offset[other_force_index]))
#define TauNew (i->tag.pointByOffset(m_omega_force_offset[force_index]))
//
//
// 4. the 3 components of the time derivative of the angular velocity
//
#define Oxdot (Tbff.x/m_J1+(m_J2-m_J3)/m_J1*Omega.y*Omega.z)
#define Oydot (Tbff.y/m_J2+(m_J3-m_J1)/m_J2*Omega.z*Omega.x)
#define Ozdot (Tbff.z/m_J3+(m_J1-m_J2)/m_J3*Omega.x*Omega.y)
//
//
// 5. then we define \dot{q} according to equation (3) in Omelyan
//
#define Q0dot  (-0.5*( Q1*Omega.x+Q2*Omega.y+Q3*Omega.z))
#define Q1dot  ( 0.5*( Q0*Omega.x-Q3*Omega.y+Q2*Omega.z))
#define Q2dot  ( 0.5*( Q3*Omega.x+Q0*Omega.y-Q1*Omega.z))
#define Q3dot  ( 0.5*(-Q2*Omega.x+Q1*Omega.y+Q0*Omega.z))
#define QdotSQ (Q0dot*Q0dot+Q1dot*Q1dot+Q2dot*Q2dot+Q3dot*Q3dot)
// 
//
// 6. then we define \ddot{q} according to equation (4) in Omelyan
//
#define Q0ddot  (-0.5*(( Q1dot*Omega.x+Q2dot*Omega.y+Q3dot*Omega.z)+(Q1*Oxdot+Q2*Oydot+Q3*Ozdot)))
#define Q1ddot  ( 0.5*(( Q0dot*Omega.x-Q3dot*Omega.y+Q2dot*Omega.z)+(Q0*Oxdot-Q3*Oydot+Q2*Ozdot)))
#define Q2ddot  ( 0.5*(( Q3dot*Omega.x+Q0dot*Omega.y-Q1dot*Omega.z)+(Q3*Oxdot+Q0*Oydot-Q1*Ozdot)))
#define Q3ddot  ( 0.5*((-Q2dot*Omega.x+Q1dot*Omega.y+Q0dot*Omega.z)+(-Q2*Oxdot+Q1*Oydot+Q0*Ozdot)))
#define QddotSQ (Q0ddot*Q0ddot+Q1ddot*Q1ddot+Q2ddot*Q2ddot+Q3ddot*Q3ddot)
// 
//
#define QdotQddot (Q0dot*Q0ddot+Q1dot*Q1ddot+Q2dot*Q2ddot+Q3dot*Q3ddot)
//
// 8. The torque saved in the tag
//
#define Tbff (i->tag.pointByOffset(m_torque_offset))
// then the timestep
//
#define DT (m_dt)
//
//
// The Lagrangian multiplier
//
//
//  if(m_first_sweep == true){
//    FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, m_colour, this,
//      Tbff = RT*TauNew;
//      );
//    m_first_sweep = false;
//  }
  FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, m_colour, this,
      Tbff = RT*TauNew;
      LAMBDA = (1.e0-QdotSQ*DT*DT/2.e0-sqrt(1.e0-QdotSQ*DT*DT-QdotQddot*DT*DT*DT-(QddotSQ-QdotSQ*QdotSQ)*DT*DT*DT*DT/4.e0));
      Q0 = Q0 + Q0dot*DT + Q0ddot*DT*DT/2.0 - LAMBDA*Q0;
      Q1 = Q1 + Q1dot*DT + Q1ddot*DT*DT/2.0 - LAMBDA*Q1;
      Q2 = Q2 + Q2dot*DT + Q2ddot*DT*DT/2.0 - LAMBDA*Q2;
      Q3 = Q3 + Q3dot*DT + Q3ddot*DT*DT/2.0 - LAMBDA*Q3;
//      if(i->mySlot == 13) MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "NormQ " << sqrt(Q0*Q0+Q1*Q1+Q2*Q2+Q3*Q3));
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
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "Colour = " << c);
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "[q0,q1,q2,q3] = [" << Q0 << ", " << Q1 << ", " << Q2 << ", " << Q3 <<"]");
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "omega = " << Omega);
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "rotmat  = " << R);
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "rotmatT = " << RT);
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "Tbff = " << Tbff);
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "TauOld = " << TauOld);
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep1()", "TauNew = " << TauNew);
      );
}


void IntegratorOmelyanNR::integrateStep2(){
  Phase *phase = M_PHASE;
  double X, DX; 
  double Y, DY; 
  double Z, DZ;
  double IX, IY, IZ;
  double FX, FY, FZ;
  double A, B, C, D, F, G, DETF;
  double GXX, GXY, GXZ, GYX, GYY, GYZ, GZX, GZY, GZZ;
  double PFX, PFY, PFZ;
  double relaerr = 1.0, residue = 1.0, norm=1.0;
  int iter, imax=1;
  point_t TauOldPrime, TauNewPrime;
//
//
// This is to define the index to the new force which is force_index 
// and the index to the old one is given as the modulo function below.
// General form for arbitrary FORCE_HIST_SIZE should read: 
// (force_index+FORCE_HIST_SIZE-1)%(FORCE_HIST_SIZE) NOTE the % operator!!!
// This only holds for calculating the previous index.
// This is sufficient for calculating all history iteratively. I you
// only need every nth historic value -> do it on your own.
// Discussion David&Andreas 14.6.2013
//
  size_t force_index = M_CONTROLLER->forceIndex();
  size_t other_force_index = (force_index+1)&(FORCE_HIST_SIZE-1);
//  MSG_DEBUG("IntegratorOmelyanNR::integrateStep2()", "force_index: " << force_index << ", other_force_index: " << other_force_index);       
//
//MSG_DEBUG("IntegratorOmelyanNR::integrateStep2()", "other_force_index " << other_force_index << " force_index " << force_index);
//exit(0);
  FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, m_colour, this,
    PFX = DT*(m_J2-m_J3)/m_J1;// /2.0;
    PFY = DT*(m_J3-m_J1)/m_J2;// /2.0;
    PFZ = DT*(m_J1-m_J2)/m_J3;// /2.0;
    TauNewPrime=RT*TauNew;
    TauOldPrime=RT*TauOld;
    IX = Omega.x+(TauOldPrime.x+TauNewPrime.x)*DT/m_J1/2.0;//+PFX*Omega.y*Omega.z;
    IY = Omega.y+(TauOldPrime.y+TauNewPrime.y)*DT/m_J2/2.0;//+PFY*Omega.z*Omega.x;
    IZ = Omega.z+(TauOldPrime.z+TauNewPrime.z)*DT/m_J3/2.0;//+PFZ*Omega.x*Omega.y;
    X = Omega.x;
    Y = Omega.y;
    Z = Omega.z;
//
// We have to solve the system of equations of second order that looks like rhis
//
//       OxNew = OxOld + (TauxOld+TauxNew)*DT/m_J1/2.0 + (m_J2-m_J3)/m_J1 * (Omega.y*Omega.z+OyNew*OzNew);
//       OyNew = OyOld + (TauyOld+TauyNew)*DT/m_J2/2.0 + (m_J3-m_J1)/m_J2 * (Omega.z*Omega.x+OzNew*OxNew);
//       OzNew = OzOld + (TauzOld+TauzNew)*DT/m_J3/2.0 + (m_J1-m_J2)/m_J3 * (Omega.x*Omega.y+OxNew*OyNew);
//
    iter=0;
    do{
     iter+=1;
//
     FX = (X-PFX*Y*Z-IX);
     FY = (Y-PFY*Z*X-IY);
     FZ = (Z-PFZ*X*Y-IZ);
     A = -PFX*Z;
     B = -PFX*Y;
     C = -PFY*Z;
     D = -PFY*X;
     F = -PFZ*Y;
     G = -PFZ*X;
     DETF = 1.0-A*C-B*F-D*G+A*D*F+B*C*G;
     GXX = (1.0-D*G)/DETF;
     GXY = (B*G-A)/DETF;
     GXZ = (A*D-B)/DETF;
     GYX = (D*F-C)/DETF;
     GYY = (1.0-B*F)/DETF;
     GYZ = (B*C-D)/DETF;
     GZX = (C*G-F)/DETF;
     GZY = (A*F-G)/DETF;
     GZZ = (1.0-A*C)/DETF;
//
     DX = -GXX*FX-GXY*FY-GXZ*FZ;
     DY = -GYX*FX-GYY*FY-GYZ*FZ;
     DZ = -GZX*FX-GZY*FY-GZZ*FZ;
//
     X = X+DX;
     Y = Y+DY;
     Z = Z+DZ;
//
//
     residue = sqrt(FX*FX+FY*FY+FZ*FZ);
     relaerr = sqrt(DX*DX+DY*DY+DZ*DZ)/sqrt(X*X+Y*Y+Z*Z);
     if(iter>10){
       MSG_DEBUG("IntegratorOmelyanNR::integrateStep2()", "iter = " << iter << ", residue = " << residue << ", relaerr = " << relaerr);
       exit(0);
     }
/*     MSG_DEBUG("IntegratorOmelyanNR::integrateStep(2)", "TauOld.x = " << TauOld.x << ", " << "TauOld.y = " << TauOld.y << ", " << "TauOld.z = " << TauOld.z << endl
	 << "TauNew.x = " << TauNew.x << ", " << "TauyNew = " << TauyNew << ", " << "TauzNew = " << TauzNew << endl);
*/
   }
   while(residue > 1.e-10 || relaerr > 1.e-10);
//  MSG_DEBUG("IntegratorOmelyanNR::integrateStep(2)", "exit do-while loop");      
//
//
   Omega.x = X;
   Omega.y = Y;
   Omega.z = Z;
   Tbff = TauNewPrime;
   if(m_normalize){
     norm = sqrt(Q0*Q0+Q1*Q1+Q2*Q2+Q3*Q3);
     Q0 = Q0/norm;
     Q1 = Q1/norm;
     Q2 = Q2/norm;
     Q3 = Q3/norm;
   }
   if(iter > imax) imax=iter;
      );
}

#ifdef _OPENMP
string IntegratorOmelyanNR::dofIntegr() {
  return m_omega_name;
}

void IntegratorOmelyanNR::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    for (int i = 0; i < SPACE_DIMS; ++i) {
      p->tag.pointByOffset(m_omega_force_offset[force_index])[i] += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i];
      (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos + i] = 0;
    }
  }
}

#endif
