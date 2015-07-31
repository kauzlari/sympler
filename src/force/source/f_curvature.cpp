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


#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "quintet.h"
#include "f_curvature.h"

// #include "valgrind/memcheck.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

#define PI M_PI

const GenFTypeConcr<FCurvature> curvature_force("FCurvature");

//---- Constructors/Destructor ----

FCurvature::FCurvature()
{
	throw gError("FCurvature::FCurvature(default)", "Should not be called. Contact the programmer.");
}

FCurvature::FCurvature(Simulation *simulation): GenQuintet(simulation)
{
	init();
}


FCurvature::~FCurvature()
{
}


void FCurvature::init()
{
	m_properties.setClassName("FCurverture");

	m_properties.setDescription(
			"This is a Force based on the curverture of the surface represented by five particles. The five particle are assumed to form a Surface that can be expressed by LangrangePolymers. This potential is only valid in R3.  F = Grad( Integral (2H^2) dA).");

	DOUBLEPC(K, m_kappa, -1,"Bending rigidity (Energy).");
	DOUBLEPC(Awork, m_AWork, -1,"Working Area (L^2).");
	
	m_kappa = 0;
	m_AWork = 0;

  m_is_pair_force = false;
  m_is_particle_force = false;
}


//---- Methods ----

void FCurvature::computeForces(int force_index)
{
  //   MSG_DEBUG("FAngular::computeForces", "START");
  // FIXME: not so nice to do that here, but currently the only possibility. See also GenQuintetr::setupAfterParticleCreation(). Is the problem solved if this function will be called for each triplet (when parallelised)?  
  if(!m_QuintetList) {
    MSG_DEBUG("FCurvature::computeForces", "requesting quintet list");
    m_QuintetList = M_PHASE -> returnQuintetList(m_force_name);
    MSG_DEBUG("FCurvature::computeForces", "got quintet list");
}

  point_t boxSize = M_PHASE->boundary()->boundingBox().size();
//   double size;

// Old: 1 line
//   for (quintetListItr p = m_QuintetList.begin(); p != m_QuintetList.end(); ++p) 
  for (quintetListItr p = m_QuintetList->begin(); p != m_QuintetList->end(); ++p) 
    {		
      point_t b00, b20, b02, b22, r11; // position vectors
      point_t nv;

      point_t GradE00, GradE20, GradE02, GradE22;
      point_t GradF00, GradF20, GradF02, GradF22;
      point_t GradG00, GradG20, GradG02, GradG22;
      point_t Ss, St; 	

      point_t GradN00, GradN20, GradN02, GradN22, GradN11;
      point_t GradM00, GradM20, GradM02, GradM22;
      point_t GradA00, GradA20, GradA02, GradA22;

      point_t F00, F20, F02, F22, F11, Null; // ForceVectors   

      double E = 0;
      double F = 0;
      double G = 0;

      double N = 0;
      double M = 0;
      double L = 0;

      int _1 = -1;
      int _j = 0;
      int _k = 0;
 
      point_t force_p00, force_p20,force_p22, force_p02, force_p11,force_p11_1;

      for (int _i = 0; _i < SPACE_DIMS; _i++)
	{   
	  b00[_i] = p->p00->r[_i] - p->p11->r[_i];
	  b20[_i] = p->p20->r[_i] - p->p11->r[_i];
	  b02[_i] = p->p02->r[_i] - p->p11->r[_i];
	  b22[_i] = p->p22->r[_i] - p->p11->r[_i];
	  r11[_i] = p->p11->r[_i];

	  // periodic BCs
	  if(b00[_i] > 0.5*boxSize[_i]) b00[_i] -= boxSize[_i]; 
	  if(b00[_i] < -0.5*boxSize[_i]) b00[_i] += boxSize[_i]; 
	  if(b20[_i] > 0.5*boxSize[_i]) b20[_i] -= boxSize[_i]; 
	  if(b20[_i] < -0.5*boxSize[_i]) b20[_i] += boxSize[_i];
 	  if(b02[_i] > 0.5*boxSize[_i]) b02[_i] -= boxSize[_i]; 
	  if(b02[_i] < -0.5*boxSize[_i]) b02[_i] += boxSize[_i]; 
	  if(b22[_i] > 0.5*boxSize[_i]) b22[_i] -= boxSize[_i]; 
	  if(b22[_i] < -0.5*boxSize[_i]) b22[_i] += boxSize[_i]; 

	  //Ss and St
	  Ss[_i] = (b22[_i]+b20[_i]-b02[_i]-b00[_i])/4;
	  St[_i] = (b22[_i]-b20[_i]+b02[_i]-b00[_i])/4;

	  // Scalars
	  E += Ss[_i]*Ss[_i]; //Ss.Ss
	  F += Ss[_i]*St[_i]; //Ss.St
	  G += St[_i]*St[_i]; //St.St

	  GradE00[_i] = (b00[_i]+b02[_i]-b20[_i]-b22[_i])/8;
	  GradE20[_i] = (b22[_i]+b20[_i]-b02[_i]-b00[_i])/8;
	  GradE02[_i] = (b00[_i]+b02[_i]-b20[_i]-b22[_i])/8; //= GradE00
	  GradE22[_i] = (b22[_i]+b20[_i]-b02[_i]-b00[_i])/8; //= GradE20
	  
	  GradF00[_i] = (b00[_i]-b22[_i])/4;
	  GradF20[_i] = (b02[_i]-b20[_i])/4;
	  GradF02[_i] = (b20[_i]-b02[_i])/4;
	  GradF22[_i] = (b22[_i]-b00[_i])/4;

	  GradG00[_i] = (b00[_i]-b02[_i]+b20[_i]-b22[_i])/8;
	  GradG20[_i] = (b00[_i]-b02[_i]+b20[_i]-b22[_i])/8; //= GradG00
	  GradG02[_i] = (b22[_i]-b20[_i]+b02[_i]-b00[_i])/8;
	  GradG22[_i] = (b22[_i]-b20[_i]+b02[_i]-b00[_i])/8; //= GradG02

	}

      // Loop for the quinatities that contain crossproducts
      for (int _i = 0; _i < SPACE_DIMS; _i++)
	{   
	  //Include crossproduct Terms		
	  _j = _i + 1;
	  _k = _i + 2;	  
 	  _1 = _1*-1;
	  if (_j > SPACE_DIMS-1) _j -= SPACE_DIMS;
	  if (_k > SPACE_DIMS-1) _k -= SPACE_DIMS;

	  //calculation of the Normalvector
	  nv[_i] = _1*(b20[_j]*(b22[_k]-b00[_k])+b02[_j]*(b00[_k]-b22[_k])-(b00[_j]-b22[_j])*(b02[_k]-b20[_k]));

	  // GradNij == GradLij
	  GradN00[_i] = _1*(b22[_j]*(b02[_k]-b20[_k])+b20[_j]*(b02[_k]+b22[_k])-b02[_j]*(b20[_k] + b22[_k]))/8;
	  GradN20[_i] = _1*(b22[_j]*(b00[_k]+b02[_k])+b02[_j]*(b00[_k]-b22[_k])-b00[_j]*(b02[_k] + b22[_k]))/8;
	  GradN02[_i] = _1*(b00[_j]*(b20[_k]+b22[_k])+b20[_j]*(b22[_k]-b00[_k])-b22[_j]*(b00[_k] + b20[_k]))/8;
	  GradN22[_i] = _1*(b00[_j]*(b20[_k]-b02[_k])+b02[_j]*(b00[_k]+b20[_k])-b20[_j]*(b00[_k] + b02[_k]))/8;
	  // Controll Force
	  GradN11[_i] = _1*((b00[_j]-b22[_j])*(b02[_k]-b20[_k])+b20[_j]*(b00[_k]-b22[_k])+b02[_j]*(b22[_k]-b00[_k]))/4;

	  GradM00[_i] = _1*(b22[_j]*(b02[_k]-b20[_k])+b02[_j]*(b20[_k]-b22[_k])+b20[_j]*(b22[_k]-b02[_k]))/16;
	  GradM20[_i] = _1*(b22[_j]*(b00[_k]-b02[_k])+b00[_j]*(b02[_k]-b22[_k])+b02[_j]*(b22[_k]-b00[_k]))/16;
	  GradM02[_i] = _1*(b22[_j]*(b20[_k]-b00[_k])+b20[_j]*(b00[_k]-b22[_k])+b00[_j]*(b22[_k]-b20[_k]))/16;
	  GradM22[_i] = _1*(b20[_j]*(b02[_k]-b00[_k])+b02[_j]*(b00[_k]-b20[_k])+b00[_j]*(b20[_k]-b02[_k]))/16;

	  //Scalars
	  N += (b22[_i]+b20[_i]+b02[_i]+b00[_i])*nv[_i]/16; //Sss.N = Sss.(Ss x St), N == L
	  M += (b22[_i]-b20[_i]-b02[_i]+b00[_i])*nv[_i]/32;//Sst.N = Sst.(Ss x St)

	  //L += (b22[_i]+b20[_i]+b02[_i]+b00[_i])*nv[_i]/16//Stt.N = Stt.(Ss x St)
	}	

	//MSG_DEBUG("FCurvature::calculate force", "calc b00" <<  b00);
	//MSG_DEBUG("FCurvature::calculate force", "calc b20" <<  b20);
	//MSG_DEBUG("FCurvature::calculate force", "calc b22" <<  b22);
	//MSG_DEBUG("FCurvature::calculate force", "calc b02" <<  b02);
	
	L = N;

	double A = E*N + G*L - 2*F*M;
	double B = E*G-F*F; 
	double H2 = A*A/B/B;

//	MSG_DEBUG("FCurvature::calculate force", "Normal Vector " <<  nv);
//	MSG_DEBUG("FCurvature::calculate force", "calc A " <<  A);
//	MSG_DEBUG("FCurvature::calculate force", "calc B " <<  B);
//	MSG_DEBUG("FCurvature::calculate force", "calc H2 " <<  H2);

	GradA00 = GradE00*N+E*GradN00+GradG00*L+G*GradN00-2*M*GradF00-2*F*GradM00;
	GradA02 = GradE02*N+E*GradN02+GradG02*L+G*GradN02-2*M*GradF02-2*F*GradM02;
	GradA22 = GradE22*N+E*GradN22+GradG22*L+G*GradN22-2*M*GradF22-2*F*GradM22;
	GradA20 = GradE20*N+E*GradN20+GradG20*L+G*GradN20-2*M*GradF20-2*F*GradM20;

	if(B == 0)
	 {
	   F00 = Null;
	   F02 = Null;
	   F22 = Null;
	   F20 = Null;	
	   F11 = Null; 
	}
	else
	 {
	   F00 = 3*(A*A)/(B*B*B*B)*(GradE00*G+E*GradG00-2*F*GradF00) - 2*A/(B*B*B)*(GradA00);
	   F02 = 3*(A*A)/(B*B*B*B)*(GradE02*G+E*GradG02-2*F*GradF02) - 2*A/(B*B*B)*(GradA02);
	   F22 = 3*(A*A)/(B*B*B*B)*(GradE22*G+E*GradG22-2*F*GradF22) - 2*A/(B*B*B)*(GradA22);
	   F20 = 3*(A*A)/(B*B*B*B)*(GradE20*G+E*GradG20-2*F*GradF20) - 2*A/(B*B*B)*(GradA20);   
	   F11 = -2*A/(B*B*B)*(GradN11*(E+G));
	}

//       MSG_DEBUG("FCurvature::calculate force", "particle " << p->b->mySlot << " cos(theta) = " << cos_a);
//       MSG_DEBUG("FCurvature::calculate force", "particle " << p->b->mySlot << " m_thetaEq = " << m_thetaEq);
//       MSG_DEBUG("FCurvature::calculate force", "particle " << p->b->mySlot << " cosEQ = " << cos(PI*m_thetaEq/180));
//       MSG_DEBUG("FCurvature::calculate force", "particle " << p->b->mySlot << " m_cosEq = " << m_cosEq);
      
      
      force_p00 = m_kappa*m_AWork/4*F00;
        //MSG_DEBUG("FCurvature::calculate force", "particle " << p->p00->mySlot << " force_p00 " << force_p00);
      force_p20 = m_kappa*m_AWork/4*F20;
        //MSG_DEBUG("FCurvature::calculate force", "particle " << p->p20->mySlot << " force_p20 " << force_p20);
      force_p22 = m_kappa*m_AWork/4*F22;
       //MSG_DEBUG("FCurvature::calculate force", "particle " << p->p22->mySlot << " force_p22 " << force_p22);
      force_p02 = m_kappa*m_AWork/4*F02;
        //MSG_DEBUG("FCurvature::calculate force", "particle " << p->p02->mySlot << " force_p02 " << force_p02);   
      force_p11   = -1*(force_p00+force_p20+force_p02+force_p22) ;
      force_p11_1 = m_kappa*m_AWork/4*F11;
       //MSG_DEBUG("FCurvature::calculate force", "particle " << p->p11->mySlot << " force_p11 " << force_p11);
       //MSG_DEBUG("FCurvature::calculate force", "particle " << p->p11->mySlot << " force_p11_1 " << force_p11_1);
      p->p00->force[force_index] += force_p00;
      p->p02->force[force_index] += force_p20;
      p->p22->force[force_index] += force_p22;
      p->p02->force[force_index] += force_p02;
      p->p11->force[force_index] += force_p11;
    }
}


#ifndef _OPENMP
void FCurvature::computeForces(Pairdist* pair, int force_index)
#else
void FCurvature::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("FCurvature::computeForces", "Fatal error: do not call FCurvature::computeForces(Pairdist* pair, int force_index)!!! Please contact the programmer!");
}


void FCurvature::computeForces(Particle* part, int force_index)
{
  throw gError("FCurvature::computeForces", "Fatal error: do not call FCurvature::computeForces(Particle* part, int force_index)!!! Please contact the programmer!");
}


void FCurvature::setup()
{
  GenQuintet::setup();
}

