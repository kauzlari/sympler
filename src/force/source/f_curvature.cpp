/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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
	m_properties.setClassName("FCurvature");

	m_properties.setDescription(
			"This is a force based on the curvature of the surface represented by five particles. The five particle are assumed to form a surface that can be expressed by lagrange polynomials. This potential is only valid in R3. Fij = kA/4*Grad_rij(H^2).");

	DOUBLEPC(k, m_k, -1,"Bending rigidity (Energy).");
	
	DOUBLEPC(C0, C0, -1,"Spontaneous curvature (1/L).");

	BOOLPC(periodic, m_periodic, "Should periodic boundary conditions be applied to the connections?");

	/*
	 * m_kappa_A Effective Bending rigidity Energy area product
	*/
	m_k = 0;

	/*
	 * C0 Spontaneous curvature used for Bilayer Simulation
	*/
	C0 = 0;

	/*!
	 * m_periodic use periodic boundarys. (03.08.2015) Cannot distiguish directions yet (implementation)
	*/
	m_periodic = true; 

	/*!
	 * This Force is neither a particle force nor a pair force!
	*/
	m_is_pair_force = false;
  	m_is_particle_force = false;
}


//---- Methods ----

void FCurvature::computeForces(int force_index)
{
  // This Function is only valis in R3
  if(SPACE_DIMS < 3) 
    {
      throw gError("FCurvature::computeForces", "This Force is only valid in R3 ");
  }


  //   MSG_DEBUG("FAngular::computeForces", "START");
  // FIXME: not so nice to do that here, but currently the only possibility. See also GenQuintetr::setupAfterParticleCreation(). Is the problem solved if this function will be called for each triplet (when parallelised)?  
  if(!m_QuintetList) {
    MSG_DEBUG("FCurvature::computeForces", "requesting quintet list");
    m_QuintetList = M_PHASE -> returnQuintetList(m_force_name);
    MSG_DEBUG("FCurvature::computeForces", "got quintet list");
  }

  point_t boxSize = M_PHASE->boundary()->boundingBox().size();

  for (quintetListItr q = m_QuintetList->begin(); q != m_QuintetList->end(); ++q) 
    {		
      point_t b00, b20, b02, b22; // position vectors
      point_t nv;

      point_t GradE00, GradE22, GradE20, GradE02;
      point_t GradF00, GradF20, GradF02, GradF22;
      point_t GradG00, GradG22, GradG20, GradG02;
      point_t Ss, St; 	

      point_t GradN00, GradN20, GradN02, GradN22,GradN11; //,GradN11;
      point_t GradM00, GradM20, GradM02, GradM22;
      point_t GradA00, GradA20, GradA02, GradA22;

      point_t GradAw1_00, GradAw4_00, GradAw1_02, GradAw2_02, GradAw2_22, GradAw3_22, GradAw3_20, GradAw4_20;
      point_t GradAw00, GradAw20, GradAw02, GradAw22; //Gradient in Area 

      point_t F00, F20, F02, F22,F11;// ForceVectors   
      
      double b00b00 = 0;
      double b20b20 = 0;
      double b02b02 = 0; 
      double b22b22 = 0; //Scalar Products
      double b00b20 = 0;
      double b20b22 = 0;
      double b22b02 = 0;
      double b02b00 = 0;

      double E = 0;
      double F = 0;
      double G = 0;

      double N = 0;
      double M = 0;
      double L = 0;

      //int _1 = -1;
      int _j = 0;
      int _k = 0;
 
      point_t force_p00, force_p20,force_p22, force_p02, force_p11,force_p11_1;

      for (int _i = 0; _i < SPACE_DIMS; _i++)
	{   
	  b00[_i] = q->p00->r[_i] - q->p11->r[_i];
	  b20[_i] = q->p20->r[_i] - q->p11->r[_i];
	  b02[_i] = q->p02->r[_i] - q->p11->r[_i];
	  b22[_i] = q->p22->r[_i] - q->p11->r[_i];

	  // periodic BCs
	  if(m_periodic == true)
	    {
	      if(b00[_i] > 0.5*boxSize[_i]) b00[_i] -= boxSize[_i]; 
	      if(b00[_i] < -0.5*boxSize[_i]) b00[_i] += boxSize[_i]; 
	      if(b20[_i] > 0.5*boxSize[_i]) b20[_i] -= boxSize[_i]; 
	      if(b20[_i] < -0.5*boxSize[_i]) b20[_i] += boxSize[_i];
 	      if(b02[_i] > 0.5*boxSize[_i]) b02[_i] -= boxSize[_i]; 
	      if(b02[_i] < -0.5*boxSize[_i]) b02[_i] += boxSize[_i]; 
	      if(b22[_i] > 0.5*boxSize[_i]) b22[_i] -= boxSize[_i]; 
	      if(b22[_i] < -0.5*boxSize[_i]) b22[_i] += boxSize[_i]; 
	  }

	  /*!
	   * Calculation of the scalar Products for the Awork calculation
	  */
	  b00b00+=b00[_i]*b00[_i];
	  b20b20+=b20[_i]*b20[_i];
	  b02b02+=b02[_i]*b02[_i];
	  b22b22+=b22[_i]*b22[_i];

	  b00b20+=b00[_i]*b20[_i];
	  b20b22+=b20[_i]*b22[_i];
	  b22b02+=b22[_i]*b02[_i];
	  b02b00+=b02[_i]*b00[_i];

	  /*!
	   * Calculation of the derivatives dS/ds and dS/dt of the Surface S
	  */
	  Ss[_i] = (b22[_i]+b20[_i]-b02[_i]-b00[_i])/4.0;
	  St[_i] = (b22[_i]-b20[_i]+b02[_i]-b00[_i])/4.0;

	  /*!
	   * The Coefficients of the first fundamental form
	  */
	  E += Ss[_i]*Ss[_i]; //Ss.Ss
	  F += Ss[_i]*St[_i]; //Ss.St
	  G += St[_i]*St[_i]; //St.St
 
	  /*!
	   *  The Gradients of the first fundamental form coefficients with respect to the paricle ij
	  */
	  GradE00[_i] = (b00[_i]+b02[_i]-b20[_i]-b22[_i])/8.0;
	  GradE20[_i] = (b22[_i]+b20[_i]-b02[_i]-b00[_i])/8.0;  //= GradE22
	  GradE02[_i] = (b00[_i]+b02[_i]-b20[_i]-b22[_i])/8.0; //= GradE00
	  GradE22[_i] = (b22[_i]+b20[_i]-b02[_i]-b00[_i])/8.0; 
	  
	  GradF00[_i] = (b00[_i]-b22[_i])/8.0;
	  GradF20[_i] = (b02[_i]-b20[_i])/8.0;
	  GradF02[_i] = (b20[_i]-b02[_i])/8.0;
	  GradF22[_i] = (b22[_i]-b00[_i])/8.0;

	  GradG00[_i] = (b00[_i]-b02[_i]+b20[_i]-b22[_i])/8.0;
	  GradG20[_i] = (b00[_i]-b02[_i]+b20[_i]-b22[_i])/8.0; //= GradG00
	  GradG02[_i] = (b22[_i]-b20[_i]+b02[_i]-b00[_i])/8.0;  //= GradG22
	  GradG22[_i] = (b22[_i]-b20[_i]+b02[_i]-b00[_i])/8.0; 

	}

  /*!
   *   Loop for alle the quantities that contain crossproducts 
  */

      for (int _i = 0; _i < SPACE_DIMS; _i++)
	{   		
	  /*!
	   *Include crossproduct Terms
	  */
	  _j = _i + 1;
	  _k = _i + 2;	  
 	  //_1 = _1*-1;
	  if (_j > SPACE_DIMS-1) _j -= SPACE_DIMS;
	  if (_k > SPACE_DIMS-1) _k -= SPACE_DIMS;

	  /*!
	   *Calculation of the non normed normalvector of the surface N = (St x Ss)
	  */
	  nv[_i]=(b02[_j]-b20[_j])*(b00[_k]-b22[_k])-(b00[_j]-b22[_j])*(b02[_k]-b20[_k]);

	  /*!
	   * The Coefficients of the second fundamental form: where N == L
	  */
	  N += (b22[_i]+b20[_i]+b02[_i]+b00[_i])*nv[_i]/16.0; //Sss.N = Sss.(Ss x St), N == L (Durch 2 ???)
	  M += (b22[_i]-b20[_i]-b02[_i]+b00[_i])*nv[_i]/32.0;//Sst.N = Sst.(Ss x St)
	  //L += (b22[_i]+b20[_i]+b02[_i]+b00[_i])*nv[_i]/16//Stt.N = Stt.(Ss x St)

	  /*!
	   * The Gradients of the second fundamental form coefficients with respect to the paricle ij
	  */
	  // GradNij == GradLij
	  GradN00[_i] = (b22[_j]*(b02[_k]-b20[_k])+b20[_j]*(b02[_k]+b22[_k])-b02[_j]*(b20[_k]+b22[_k]))/8.0;
	  GradN20[_i] = (b22[_j]*(b00[_k]+b02[_k])+b02[_j]*(b00[_k]-b22[_k])-b00[_j]*(b02[_k]+b22[_k]))/8.0;
	  GradN02[_i] = (b00[_j]*(b20[_k]+b22[_k])+b20[_j]*(b22[_k]-b00[_k])-b22[_j]*(b00[_k]+b20[_k]))/8.0;
	  GradN22[_i] = (b00[_j]*(b20[_k]-b02[_k])+b02[_j]*(b00[_k]+b20[_k])-b20[_j]*(b00[_k]+b02[_k]))/8.0;
	  // Controll Force
	  //GradN11[_i] = ((b00[_j]-b22[_j])*(b02[_k]-b20[_k])+b20[_j]*(b00[_k]-b22[_k])+b02[_j]*(b22[_k]-b00[_k]))/4.0;

	  GradM00[_i] = (b22[_j]*(b02[_k]-b20[_k])+b02[_j]*(b20[_k]-b22[_k])+b20[_j]*(b22[_k]-b02[_k]))/16.0;
	  GradM20[_i] = (b22[_j]*(b00[_k]-b02[_k])+b00[_j]*(b02[_k]-b22[_k])+b02[_j]*(b22[_k]-b00[_k]))/16.0;
	  GradM02[_i] = (b22[_j]*(b20[_k]-b00[_k])+b20[_j]*(b00[_k]-b22[_k])+b00[_j]*(b22[_k]-b20[_k]))/16.0;
	  GradM22[_i] = (b20[_j]*(b02[_k]-b00[_k])+b02[_j]*(b00[_k]-b20[_k])+b00[_j]*(b20[_k]-b02[_k]))/16.0;

	  /*!
	   * Gradients of the Area with respect to the particle ij , 2Grad[Awork]=2GradAw/(Awork*4)/8=2GradAw_1/Awork1
	   */	  	
	  GradAw1_00[_i] = (b00[_i]*(b02[_j]*b02[_j]+b02[_k]*b02[_k])-b02[_i]*(b00[_j]*b02[_j]+b00[_k]*b02[_k]));
	  GradAw4_00[_i] = (b00[_i]*(b20[_j]*b20[_j]+b20[_k]*b20[_k])-b20[_i]*(b00[_j]*b20[_j]+b00[_k]*b20[_k]));

	  GradAw1_02[_i] = (b02[_i]*(b00[_j]*b00[_j]+b00[_k]*b00[_k])-b00[_i]*(b02[_j]*b00[_j]+b02[_k]*b00[_k]));
	  GradAw2_02[_i] = (b02[_i]*(b22[_j]*b22[_j]+b22[_k]*b22[_k])-b22[_i]*(b02[_j]*b22[_j]+b02[_k]*b22[_k]));

	  GradAw2_22[_i] = (b22[_i]*(b02[_j]*b02[_j]+b02[_k]*b02[_k])-b02[_i]*(b22[_j]*b02[_j]+b22[_k]*b02[_k]));
	  GradAw3_22[_i] = (b22[_i]*(b20[_j]*b20[_j]+b20[_k]*b20[_k])-b20[_i]*(b22[_j]*b20[_j]+b22[_k]*b20[_k]));

	  GradAw3_20[_i] = (b20[_i]*(b22[_j]*b22[_j]+b22[_k]*b22[_k])-b22[_i]*(b20[_j]*b22[_j]+b20[_k]*b22[_k]));
	  GradAw4_20[_i] = (b20[_i]*(b00[_j]*b00[_j]+b00[_k]*b00[_k])-b00[_i]*(b20[_j]*b00[_j]+b20[_k]*b00[_k]));

	}	
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot << " Nv" << nv);

	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot  <<  "calc GradAw1_00" <<   GradAw1_00);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot  <<  "calc GradAw4_00" <<   GradAw4_00);

	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw1_02" <<   GradAw1_02);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw2_01" <<   GradAw2_02);

	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw2_22" <<   GradAw2_22);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw3_22" <<   GradAw3_22);
	
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw3_20" <<   GradAw3_20);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw4_20" <<   GradAw4_20);
	
	

	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot  <<  " calc GradN20" <<  GradN20);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot  <<  " calc GradN02" <<  GradN02);

	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot  <<  " calc GradM20" <<  GradM20);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot  <<  " calc GradM02" <<  GradM02);

	L = N;

	// The mean curvature caculates as H = A/B*1/2 
	double A = E*N + G*L - 2*F*M;
	double B = E*G-F*F; 

	//Working Area calculation 
	double Awork1 = sqrt(b02b02*b00b00-(b02b00*b02b00));
	double Awork2 = sqrt(b22b22*b02b02-(b22b02*b22b02));
	double Awork3 = sqrt(b20b20*b22b22-(b20b22*b20b22));
	double Awork4 = sqrt(b00b00*b20b20-(b00b20*b00b20));
        double Awork =  Awork1 + Awork2 + Awork3 + Awork4;
	Awork = Awork/2;
        //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot << " Awork " << Awork);
	
	//Gradient in Working Area
	GradAw00 = (GradAw1_00/Awork1+GradAw4_00/Awork4)/2.0;
	GradAw02 = (GradAw1_02/Awork1+GradAw2_02/Awork2)/2.0;
	GradAw22 = (GradAw2_22/Awork2+GradAw3_22/Awork3)/2.0;
	GradAw20 = (GradAw3_20/Awork3+GradAw4_20/Awork4)/2.0;

	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw00" <<   GradAw00);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw02" <<   GradAw02);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw22" <<   GradAw22);
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot <<  "calc GradAw20" <<   GradAw20);

	// Variables for faster Caclulation
	double B4 = B*B*B*B;
	double B3 = B*B*B;
	double A2 = A*A;

	//Spontaeous Curvature Prefactor (H-C0/2)^2 (C0 = 0)
	double HC02 = A2/B3/4.0;
	
//	MSG_DEBUG("FCurvature::calculate force", "Normal Vector " <<  nv);
//	MSG_DEBUG("FCurvature::calculate force", "calc A " <<  A);
//	MSG_DEBUG("FCurvature::calculate force", "calc B " <<  B);
//	MSG_DEBUG("FCurvature::calculate force", "calc H2 " <<  H2);

	/*
	 *Calculation of the gradient of A (Grad(A)) Chain and Prodict rule on A
	*/
	GradA00 = GradE00*N+E*GradN00+GradG00*L+G*GradN00-2.0*M*GradF00-2.0*F*GradM00; 
	GradA22 = GradE22*N+E*GradN22+GradG22*L+G*GradN22-2.0*M*GradF22-2.0*F*GradM22;  
	GradA02 = GradE02*N+E*GradN02+GradG02*L+G*GradN02-2.0*M*GradF02-2.0*F*GradM02; //GradE02 = GradE00,GradG02 = GradG22
	//GradA02 = GradE00*N+E*GradN02+GradG22*L+G*GradN02-2*M*GradF02-2*F*GradM02;	
	GradA20 = GradE20*N+E*GradN20+GradG20*L+G*GradN20-2.0*M*GradF20-2.0*F*GradM20; //GradE20 = GradE22,GradG20 = GradG00
	//GradA20 = GradE22*N+E*GradN20+GradG00*L+G*GradN20-2*M*GradF20-2*F*GradM20; 

	/*
	 *Calculation of: -4*Grad(H^2) = -Grad(A^2/B^2) = 3*A^2/B^4*Grad(B) - 2*A/B^3*Grad(A)
	 *Where: Grad(B) = Grad(E)*G + Grad(G)*E - 2*F*Grad(F)
	*/
	F00 = (3.0*A2/B4*(GradE00*G+E*GradG00-2.0*F*GradF00)-2.0*A/B3*GradA00);
	F22 = (3.0*A2/B4*(GradE22*G+E*GradG22-2.0*F*GradF22)-2.0*A/B3*GradA22);
	F02 = (3.0*A2/B4*(GradE02*G+E*GradG02-2.0*F*GradF02)-2.0*A/B3*GradA02); //GradE02 = GradE00,GradG02 = GradG22
	//F02 = 3*A2/B4*(GradE00*G+E*GradG22-2*F*GradF02) - 2*A/B3*(GradA02);
	F20 = (3.0*A2/B4*(GradE20*G+E*GradG20-2.0*F*GradF20)-2.0*A/B3*GradA20); //GradE20 = GradE22,GradG20 = GradG00
	//F20 = 3*A2/B4*(GradE22*G+E*GradE22-2*F*GradF20) - 2*A/B3*(GradA20); 

	
	double B12 = sqrt(B);

	//double H = A/B12/B/2;
	//MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot << " H=" << H);

	if(C0 != 0 ){		
		double B32 = B12*B;		
		double B52 = B32*B;

		//Prefactor Grad(4*H*C0) = 4*c0*Grad(H)
		F00 = F00 + (2.0*C0*GradA00/B32-3.0*C0*A/B52*(GradE00*G+E*GradG00-2.0*F*GradF00))*2;
		F22 = F22 + (2.0*C0*GradA22/B32-3.0*C0*A/B52*(GradE22*G+E*GradG22-2.0*F*GradF22))*2;
		F02 = F02 + (2.0*C0*GradA02/B32-3.0*C0*A/B52*(GradE02*G+E*GradG02-2.0*F*GradF02))*2;
		F20 = F20 + (2.0*C0*GradA20/B32-3.0*C0*A/B52*(GradE20*G+E*GradG20-2.0*F*GradF20))*2;
		//Spontaeous Curvature Prefactor (H-C0/2)^2 (C0 != 0)
		HC02 = (A/B32-C0)*(A/B32-C0)/4.0;


	}
	
	//F11 = -2*A/(B*B*B)*(GradN11*(E+G));

//       MSG_DEBUG("FCurvature::calculate force", "particle " << q->b->mySlot << " cos(theta) = " << cos_a);
//       MSG_DEBUG("FCurvature::calculate force", "particle " << q->b->mySlot << " m_thetaEq = " << m_thetaEq);
//       MSG_DEBUG("FCurvature::calculate force", "particle " << q->b->mySlot << " cosEQ = " << cos(PI*m_thetaEq/180));
//       MSG_DEBUG("FCurvature::calculate force", "particle " << q->b->mySlot << " m_cosEq = " << m_cosEq);
      
	
      force_p00 = m_k*(Awork/2*F00 - 2.0*GradAw00*HC02); 
        //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p00->mySlot << " force_p00 " << force_p00);
      force_p20 = m_k*(Awork/2*F20 - 2.0*GradAw20*HC02);
        //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p20->mySlot << " force_p20 " << force_p20);
      force_p22 = m_k*(Awork/2*F22 - 2.0*GradAw22*HC02);
       //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p22->mySlot << " force_p22 " << force_p22);
      force_p02 = m_k*(Awork/2*F02 - 2.0*GradAw02*HC02);
        //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p02->mySlot << " force_p02 " << force_p02);   
      force_p11 = -1.0*(force_p00+force_p20+force_p02+force_p22);

      //Test Force consitency force_p11_1 === force_p11Grad
      //force_p11_1 = m_k*Awork/2*F11-force_p11;
      //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot << " force_dp11 " << force_p11_1);
      //MSG_DEBUG("FCurvature::calculate force", "particle " << q->p11->mySlot << " force_p11_1 " << force_p11_1);

      q->p00->force[force_index] += force_p00;
      q->p02->force[force_index] += force_p02;
      q->p22->force[force_index] += force_p22;
      q->p20->force[force_index] += force_p20;
      q->p11->force[force_index] += force_p11;
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

