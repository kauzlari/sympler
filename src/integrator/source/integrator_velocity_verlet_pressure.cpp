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
#include "integrator_velocity_verlet_pressure.h"
#include "cell.h"
#include "manager_cell.h"

#ifdef HAVE_JAMA_JAMA_LU_H

#include "jama/jama_lu.h"


using namespace std;

#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorVelocityVerletPressure>
		integrator_velocity_verlet_pressure("IntegratorVelocityVerletPressure");

//---- Constructors/Destructor ----

IntegratorVelocityVerletPressure::IntegratorVelocityVerletPressure(
		Controller *controller) :
	IntegratorVelocityVerlet(controller) {
	init();
}

IntegratorVelocityVerletPressure::~IntegratorVelocityVerletPressure() {
}

//---- Methods ----

void IntegratorVelocityVerletPressure::init() {
	// some modules need to know whether there is an Integrator,
	// which changes positions, that's why the following
	m_properties.setClassName("IntegratorPosition");
	m_properties.setName("IntegratorVelocityVerletPressure");

	m_properties.setDescription(
			"Integrates the position and momentum coordinates of each particle according to the Velocity-Verlet Algorithm");

	STRINGPC
	(nablaWF, m_nablaWF_symbol,
			"Symbol assigned to nabla of the weighting function")
		;

	STRINGPC
	(vector, m_vector_name,
			"Full name of the additional scalar field, usable as attribute in other modules")
		;

	STRINGPC
	(symbol, m_vector_symbol,
			"Symbol assigned to the gradient P, usable in algebraic expressions")
		;

	STRINGPC
	(velcorr, m_velcorr_name,
			"Full name of the additional scalar field, usable as attribute in other modules")
		;

	STRINGPC
	(velcorr_symbol, m_velcorr_symbol,
			"Symbol assigned to the gradient P, usable in algebraic expressions")
		;

	STRINGPC
	(rho, m_rho_name,
			"Symbol assigned to the vector pair interactions, usable in algebraic expressions, for the off-diagonal elements of the DIV matrix")
		;
	STRINGPC
	(tensor, m_tensor_symbol,
			"Name of the correction matrix.")
		;

}

void IntegratorVelocityVerletPressure::setup() {
	IntegratorVelocityVerlet::setup();

	m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species), M_MANAGER->getColour(
			m_species)/*m_species*/);
	m_cp->setNeedPairs(true);

	m_vector_offset = Particle::s_tag_format[m_colour].addAttribute(
			m_vector_name, DataFormat::POINT, true, m_vector_symbol).offset;

	m_velcorr_offset = Particle::s_tag_format[m_colour].addAttribute(
			m_velcorr_name, DataFormat::POINT, true, m_velcorr_symbol).offset;
}

void IntegratorVelocityVerletPressure::isAboutToStart() {

  IntegratorVelocityVerlet::isAboutToStart();
  
  Phase *phase = M_PHASE;
  
  m_rho_offset
    = Particle::s_tag_format[m_colour].attrByName(m_rho_name).offset;
  
  m_tensor_offset
    = Particle::s_tag_format[m_colour].offsetByName(m_tensor_symbol);
  
  numberOfParticles = M_PHASE->returnNofPartC(m_colour);
  grad = MArray2D(numberOfParticles, 3, 0);
  
}

void IntegratorVelocityVerletPressure::integrateStep1() {
	M_PHASE->invalidatePositions((IntegratorPosition*) this);

	FOR_EACH_FREE_PARTICLE_C
	  (M_PHASE,m_colour,
	   __iSLFE->tag.pointByOffset(m_vector_offset)[0]+=grad(__iSLFE->mySlot,0);
	   __iSLFE->tag.pointByOffset(m_vector_offset)[1]+=grad(__iSLFE->mySlot,1);
	   __iSLFE->tag.pointByOffset(m_vector_offset)[2]+=grad(__iSLFE->mySlot,2);
	);

}

void IntegratorVelocityVerletPressure::integrateStep2() {

	Phase *phase = M_PHASE;
	/*  m_force_index = M_CONTROLLER->forceIndex();*/
	size_t force_index = M_CONTROLLER->forceIndex();

	//   m_other_force_index = (m_force_index+1)&(FORCE_HIST_SIZE-1);
	size_t other_force_index = (/*m_*/force_index + 1) & (FORCE_HIST_SIZE - 1);

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
	 i->v += i->dt/2 * i->force[/*((IntegratorPosition*) data)->m_*/force_index]/m_mass;	 
	 );

 	phase->invalidateVelocities();

	//numberOfParticles = M_PHASE->returnNofPartC(m_colour);
	/*Construction of the LAPLACE COEFFICIENT MATRIX BEGIN ************************************************************/

	FOR_EACH_PAIR(m_cp,

			DataFormat::attribute_t Attr =
			m_cp->tagFormat().attrByName(m_nablaWF_symbol);
			m_nablaWF_offset = Attr.offset;

	)
		;
	//MSG_DEBUG("IntegratorPosition::integrateStep2", "numPart" <<numberOfParticles);
	/*Initialization of the gradient, one particle is considered frozen and not included in the calculation*/

	m_matrixD_ijx = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixD_ijy = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixD_ijz = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);

	m_matrixR_ijx = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixR_ijy = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixR_ijz = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);

	/*Initialization of the laplace matrix, one particle is considered frozen and not included in the calculation*/

	m_matrixL_ikx = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixL_iky = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixL_ikz = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);
	m_matrixL = MArray2D(numberOfParticles - 1, numberOfParticles - 1, 0);

	FOR_EACH_PAIR(m_cp,
			Particle* firstP = pair->firstPart();
			Particle* secondP = pair->secondPart();
			if(pair->firstPart()->mySlot<numberOfParticles-1 && pair->secondPart()->mySlot<numberOfParticles-1)
			{
				/*for grad x*/

				m_matrixD_ijx(firstP->mySlot,firstP->mySlot)+=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/secondP->tag.doubleByOffset(m_rho_offset))[0]);
				m_matrixD_ijx(secondP->mySlot,secondP->mySlot)-=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/firstP->tag.doubleByOffset(m_rho_offset))[0]);
				m_matrixD_ijx(firstP->mySlot,secondP->mySlot)-=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/secondP->tag.doubleByOffset(m_rho_offset))[0]);
				m_matrixD_ijx(secondP->mySlot,firstP->mySlot)+=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/firstP->tag.doubleByOffset(m_rho_offset))[0]);

				/*for grad y*/
				m_matrixD_ijy(firstP->mySlot,firstP->mySlot)+=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/secondP->tag.doubleByOffset(m_rho_offset))[1]);
				m_matrixD_ijy(secondP->mySlot,secondP->mySlot)-=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/firstP->tag.doubleByOffset(m_rho_offset))[1]);
				m_matrixD_ijy(firstP->mySlot,secondP->mySlot)-=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/secondP->tag.doubleByOffset(m_rho_offset))[1]);
				m_matrixD_ijy(secondP->mySlot,firstP->mySlot)+=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/firstP->tag.doubleByOffset(m_rho_offset))[1]);

				/*for grad y*/

				m_matrixD_ijz(firstP->mySlot,firstP->mySlot)+=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/secondP->tag.doubleByOffset(m_rho_offset))[2]);
				m_matrixD_ijz(secondP->mySlot,secondP->mySlot)-=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/firstP->tag.doubleByOffset(m_rho_offset))[2]);
				m_matrixD_ijz(firstP->mySlot,secondP->mySlot)-=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/secondP->tag.doubleByOffset(m_rho_offset))[2]);
				m_matrixD_ijz(secondP->mySlot,firstP->mySlot)+=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/firstP->tag.doubleByOffset(m_rho_offset))[2]);

			}
	)
		;

	FOR_EACH_PAIR(m_cp,
			Particle* firstP = pair->firstPart();
			Particle* secondP = pair->secondPart();
			if(pair->firstPart()->mySlot<numberOfParticles-1 && pair->secondPart()->mySlot<numberOfParticles-1)
			{
				/*for grad x*/

				m_matrixR_ijx(firstP->mySlot,firstP->mySlot)+=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[0]);
				m_matrixR_ijx(secondP->mySlot,secondP->mySlot)-=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[0]);
				m_matrixR_ijx(firstP->mySlot,secondP->mySlot)-=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[0]);
				m_matrixR_ijx(secondP->mySlot,firstP->mySlot)+=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[0]);

				/*for grad y*/
				m_matrixR_ijy(firstP->mySlot,firstP->mySlot)+=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[1]);
				m_matrixR_ijy(secondP->mySlot,secondP->mySlot)-=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[1]);
				m_matrixR_ijy(firstP->mySlot,secondP->mySlot)-=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[1]);
				m_matrixR_ijy(secondP->mySlot,firstP->mySlot)+=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[1]);

				/*for grad y*/

				m_matrixR_ijz(firstP->mySlot,firstP->mySlot)+=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[2]);
				m_matrixR_ijz(secondP->mySlot,secondP->mySlot)-=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[2]);
				m_matrixR_ijz(firstP->mySlot,secondP->mySlot)-=(((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[2]);
				m_matrixR_ijz(secondP->mySlot,firstP->mySlot)+=(((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))/(secondP->tag.doubleByOffset(m_rho_offset)*firstP->tag.doubleByOffset(m_rho_offset)))[2]);

			}
	)
		;
	//cout<<"m_matrixD_ijx"<<"\n";
	//cout<<(m_matrixD_ijx);
	//cout<<"m_matrixD_ijy"<<"\n";
	//cout<<(m_matrixD_ijy);
	//cout<<"m_matrixD_ijz"<<"\n";
	//cout<<(m_matrixD_ijz);


	/*Summation over the gradient matrix to get the coefficients of the Laplace matrix*/
	//for (int i=0;i<numberOfParticles-1;i++)
	//{
	//	for(int k=0;k<numberOfParticles-1;k++)
	//	{
	//		for(int j=0;j<numberOfParticles-1;j++)
	//		{
	//		m_matrixL_ikx(i,k)+=(m_matrixD_ijx(i,j)*m_matrixD_ijx(j,k));
	//		m_matrixL_iky(i,k)+=m_matrixD_ijy(i,j)*m_matrixD_ijy(j,k);
	//		m_matrixL_ikz(i,k)+=m_matrixD_ijz(i,j)*m_matrixD_ijz(j,k);
	//		}
	//	}
	//}
	//;
	for (int i = 0; i < numberOfParticles - 1; i++) {
		for (int k = 0; k < numberOfParticles - 1; k++) {
			for (int j = 0; j < numberOfParticles - 1; j++) {
				m_matrixL_ikx(i, k) += m_matrixR_ijx(i, j) * (m_matrixD_ijx(j,
						k) - m_matrixD_ijx(i, k));
				m_matrixL_iky(i, k) += m_matrixR_ijy(i, j) * (m_matrixD_ijy(j,
						k) - m_matrixD_ijy(i, k));
				m_matrixL_ikz(i, k) += m_matrixR_ijz(i, j) * (m_matrixD_ijy(j,
						k) - m_matrixD_ijz(i, k));
			}
		}
	};
	m_matrixL = m_matrixL_ikx + m_matrixL_iky + m_matrixL_ikz;

	//MArray2D velz(numberOfParticles-1,1,0);
	//MArray2D lapx(numberOfParticles-1,1,0);
	//MArray2D lapy(numberOfParticles,1,0);
	//MArray2D lapz(numberOfParticles,1,0);
	//cout<<"m_matrixikx";
	//cout<<(m_matrixL_ikx);
	//cout<<"m_matrixL_iky"<<"\n";
	//cout<<(m_matrixL_iky);
	//cout<<"m_matrixL_ikz"<<"\n";
	//cout<<(m_matrixL_ikz);


	//	FOR_EACH_FREE_PARTICLE_C(M_PHASE,m_colour,
	//if (i->mySlot<numberOfParticles-1)
	//
	//			velz(__iSLFE->mySlot,0)=__iSLFE->v[2];
	//
	//	)
	//		;
	//	for(int k=0;k<numberOfParticles-1;k++)
	//		{
	//			for(int j=0;j<numberOfParticles-1;j++)
	//			{
	//			lapx(k,0)+=m_matrixL_ikx(k,j)*velz(j,0);
	//			//lapy(k,0)+=m_matrixL_iky(k,j)*velz(j,0);
	//			//lapz(k,0)+=m_matrixL_ikz(k,j)*velz(j,0);
	//
	//			}
	//		}
	//lapx=matmult(m_matrixL_ikx,velz);
	//cout<<"\nLaplace x"<<lapx;
	////    lapy=matmult(m_matrixL_iky,velz);
	//   cout<<"\nLaplace y"<<lapy;
	////    lapz=matmult(m_matrixL_ikz,velz);
	//   cout<<"\nLaplace z"<<lapz;
	//the resulting matrix is the coefficient laplace matrix

	/*Construction of the LAPLACE COEFFICIENT MATRIX END ************************************************************/

	/*Construction of the DIV V VECTOR BEGIN ************************************************************/
	/*The divergence is calculated after calculating the grad of the velocity vector*/

	grad_v = MArray2D(numberOfParticles - 1, 9, 0);
	div_v = MArray2D(numberOfParticles - 1, 1, 0);
	FOR_EACH_FREE_PAIR(m_cp,
			Particle* firstP = pair->firstPart();
			Particle* secondP = pair->secondPart();
			if(firstP->mySlot<numberOfParticles-1 && secondP->mySlot<numberOfParticles-1)
			{
				/*grad Vx*/
				grad_v(firstP->mySlot,0)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*(firstP->v[0]-secondP->v[0])/secondP->tag.doubleByOffset(m_rho_offset))[0];
				grad_v(secondP->mySlot,0)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[0]-secondP->v[0]))/firstP->tag.doubleByOffset(m_rho_offset))[0];

				grad_v(firstP->mySlot,1)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[0]-secondP->v[0]))/secondP->tag.doubleByOffset(m_rho_offset))[1];
				grad_v(secondP->mySlot,1)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[0]-secondP->v[0]))/firstP->tag.doubleByOffset(m_rho_offset))[1];

				grad_v(firstP->mySlot,2)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[0]-secondP->v[0]))/secondP->tag.doubleByOffset(m_rho_offset))[2];
				grad_v(secondP->mySlot,2)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[0]-secondP->v[0]))/firstP->tag.doubleByOffset(m_rho_offset))[2];
				/*grad Vy*/
				grad_v(firstP->mySlot,3)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*(firstP->v[1]-secondP->v[1])/secondP->tag.doubleByOffset(m_rho_offset))[0];
				grad_v(secondP->mySlot,3)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[1]-secondP->v[1]))/firstP->tag.doubleByOffset(m_rho_offset))[0];

				grad_v(firstP->mySlot,4)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[1]-secondP->v[1]))/secondP->tag.doubleByOffset(m_rho_offset))[1];
				grad_v(secondP->mySlot,4)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[1]-secondP->v[1]))/firstP->tag.doubleByOffset(m_rho_offset))[1];

				grad_v(firstP->mySlot,5)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[1]-secondP->v[1]))/secondP->tag.doubleByOffset(m_rho_offset))[2];
				grad_v(secondP->mySlot,5)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[1]-secondP->v[1]))/firstP->tag.doubleByOffset(m_rho_offset))[2];
				/*grad Vz*/

				grad_v(firstP->mySlot,6)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*(firstP->v[2]-secondP->v[2])/secondP->tag.doubleByOffset(m_rho_offset))[0];
				grad_v(secondP->mySlot,6)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[2]-secondP->v[2]))/firstP->tag.doubleByOffset(m_rho_offset))[0];

				grad_v(firstP->mySlot,7)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[2]-secondP->v[2]))/secondP->tag.doubleByOffset(m_rho_offset))[1];
				grad_v(secondP->mySlot,7)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[2]-secondP->v[2]))/firstP->tag.doubleByOffset(m_rho_offset))[1];

				grad_v(firstP->mySlot,8)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[2]-secondP->v[2]))/secondP->tag.doubleByOffset(m_rho_offset))[2];
				grad_v(secondP->mySlot,8)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((firstP->v[2]-secondP->v[2]))/firstP->tag.doubleByOffset(m_rho_offset))[2];
			}
	)
		;

	FOR_EACH_FREE_PARTICLE_C
	  (M_PHASE,m_colour,
	   if(__iSLFE->mySlot<numberOfParticles-1)
	     {
	       (div_v)(__iSLFE->mySlot,0)=grad_v(__iSLFE->mySlot,0)+grad_v(__iSLFE->mySlot,4)+grad_v(__iSLFE->mySlot,8);
	     }
	   );

		//for each particles but the last,the attribute calculated form the input file for DIV V
		/*	div_v = MArray2D(numberOfParticles - 1, 1, 0);

		 FOR_EACH_FREE_PARTICLE_C(M_PHASE,m_colour,
		 if(i->mySlot<numberOfParticles-1)
		 {
		 (div_v)(i->mySlot,0)=i->tag.doubleByOffset(m_divergence_offset);

		 }
		 )
		 ;*/

		/*Divergence divided by delta t according to Fletcher*/
		div_v = div_v.scalarmult(1 / m_dt);

	//for each particle but the last,add the off-diagonals to the div vels vector
	//	FOR_EACH_PAIR(m_cp,
	//			Particle* firstP = pair->firstPart();
	//			Particle* secondP = pair->secondPart();
	//			if(pair->firstPart()->mySlot<numberOfParticles-1 && pair->secondPart()->mySlot<numberOfParticles-1)
	//			{
	//				/*MSG_DEBUG("IntegratorPosition::integrateStep2", "pair contribution " << (pair->tag.pointByOffset(m_DivPair_offset)));
	//				MSG_DEBUG("IntegratorPosition::integrateStep2", "1myslot " << (pair->firstPart()->mySlot));
	//				MSG_DEBUG("IntegratorPosition::integrateStep2", "2myslot " << (pair->secondPart()->mySlot));
	//				MSG_DEBUG("IntegratorPosition::integrateStep2", "1vels " << (pair->firstPart()->v));
	//				MSG_DEBUG("IntegratorPosition::integrateStep2", "2vels " << (pair->secondPart()->v));
	//				*/
	//
	//				(div_v)(pair->firstPart()->mySlot,0)+=(pair->firstPart()->v-pair->secondPart()->v)*(pair->tag.pointByOffset(m_DivPair_offset)/firstP->tag.doubleByOffset(m_rho_offset));
	//				(div_v)(pair->secondPart()->mySlot,0)+=(pair->firstPart()->v-pair->secondPart()->v)*(pair->tag.pointByOffset(m_DivPair_offset)/secondP->tag.doubleByOffset(m_rho_offset));
	//			}
	//	);

	//div_v=div_v.scalarmult(1/m_dt);

	//cout << "Div" << div_v;

	/*Construction of the DIV VECTOR END ************************************************************/

	/*Solve the LSE and write the result to another vector BEGIN*/

	JAMA::LU<double> lu(*m_matrixL.p_array2d);
	MArray2D* result = new MArray2D(numberOfParticles - 1, 1);
	result->p_array2d = new TNT::Array2D<double>(lu.solve(*(div_v.p_array2d)));
	MArray2D* result_c = new MArray2D(numberOfParticles, 1, 0);
	for (int i = 0; i < numberOfParticles - 1; i++) {
		(*result_c)(i, 0) = (*result)(i, 0);
	};
	//add the last frozen particle with laplace p of 0
	(*result_c)(numberOfParticles - 1, 0) = 0;

	/*Solve the LSE and write the result to another vector END*/

	/*Take the GRAD of the solution BEGIN*/

	FOR_EACH_FREE_PAIR(m_cp,
			Particle* firstP = pair->firstPart();
			Particle* secondP = pair->secondPart();

			grad(firstP->mySlot,0)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((*result_c)(firstP->mySlot,0)-(*result_c)(secondP->mySlot,0))/secondP->tag.doubleByOffset(m_rho_offset))[0];
			grad(secondP->mySlot,0)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((*result_c)(firstP->mySlot,0)-(*result_c)(secondP->mySlot,0))/firstP->tag.doubleByOffset(m_rho_offset))[0];
			grad(firstP->mySlot,1)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((*result_c)(firstP->mySlot,0)-(*result_c)(secondP->mySlot,0))/secondP->tag.doubleByOffset(m_rho_offset))[1];
			grad(secondP->mySlot,1)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((*result_c)(firstP->mySlot,0)-(*result_c)(secondP->mySlot,0))/firstP->tag.doubleByOffset(m_rho_offset))[1];
			grad(firstP->mySlot,2)+=((((firstP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((*result_c)(firstP->mySlot,0)-(*result_c)(secondP->mySlot,0))/secondP->tag.doubleByOffset(m_rho_offset))[2];
			grad(secondP->mySlot,2)+=((((secondP->tag.tensorByOffset(m_tensor_offset))*pair->tag.pointByOffset(m_nablaWF_offset)))*((*result_c)(firstP->mySlot,0)-(*result_c)(secondP->mySlot,0))/firstP->tag.doubleByOffset(m_rho_offset))[2];

	)
		;

	/*Take the GRAD of the solution END*/
	/*Correction of the velocity*/

	FOR_EACH_FREE_PARTICLE_C
	  (M_PHASE,m_colour,
	   
	   __iSLFE->tag.pointByOffset(m_velcorr_offset)[0]=-m_dt*grad(__iSLFE->mySlot,0)/__iSLFE->tag.doubleByOffset(m_rho_offset);
	   __iSLFE->tag.pointByOffset(m_velcorr_offset)[1]=-m_dt*grad(__iSLFE->mySlot,1)/__iSLFE->tag.doubleByOffset(m_rho_offset);
	   __iSLFE->tag.pointByOffset(m_velcorr_offset)[2]=-m_dt*grad(__iSLFE->mySlot,2)/__iSLFE->tag.doubleByOffset(m_rho_offset);
	   );

}

void IntegratorVelocityVerletPressure::integratePosition(Particle* p, Cell* cell) {
  size_t force_index;
  force_index = ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();

  point_t pt = p->force[force_index] / m_mass;
  cell->doCollision(p, p->r, p->v, pt, (IntegratorPosition*) this);
  
  p->r += p->dt * (p->v + 0.5 * p->dt * pt);
  //MSG_DEBUG("IntegratorVelocityVerlet::integratePosition", name() << "pos_force= " <<  p->force[force_index]);
  //MSG_DEBUG("IntegratorVelocityVerlet::integratePosition", name() << "position= " <<  p->r);

}

void IntegratorVelocityVerletPressure::integrateVelocity(Particle* p) {
	size_t force_index;
	force_index
			= ((Controller*) m_parent/*integrator->parent()*/)->forceIndex();
	// MSG_DEBUG("IntegratorVelocityVerlet::integrateVelocity", name() << "v BEFORE = " << p->v);
	p->v += m_lambda * (p->dt * p->force[force_index] / m_mass);
	//MSG_DEBUG("IntegratorVelocityVerlet::integrateVelocity", name() << "v AFTER = " << p->v);

}

void IntegratorVelocityVerletPressure::solveHitTimeEquation(
		WallTriangle* wallTriangle, const Particle* p, const point_t &force,
		vector<double>* results) {
	double a, b, c;
	double t0, t1;
	int n;
	point_t surface_normal = wallTriangle->normal();

	a = surface_normal * force / m_mass / 2;
	b = surface_normal * p->v;
	c = surface_normal * p->r - wallTriangle->nDotR();

	if (a != 0) {
		n = gsl_poly_solve_quadratic(a, b, c, &t0, &t1);

		//     MSG_DEBUG("IntegratorVelocityVerlet::solveHitTimeEquation", "n = " << n << ", t0 = " << t0 << ", t1 = " << t1 << " for a = " << a << ", b = " << b << ", c = " << c << ", force = " << force << ", v = " << p->v << ", r = " << p->r << ", ndotr = " << wallTriangle->nDotR() << ", surfnormal = " << surface_normal);

		if (n == 0 || (t0 < c_wt_time_eps && t1 < c_wt_time_eps)) {
		}

		else {
			if (t0 < c_wt_time_eps) {
				t0 = t1;
				results->push_back(t0);
			}

			results->push_back(t0);
			results->push_back(t1);
			sort(results->begin(), results->end());
		}
	} else {
		t0 = -c / b;

		if (t0 < c_wt_time_eps) {
			n = 0;
		}

		n = 1;
		results->push_back(t0);
	}
}

void IntegratorVelocityVerletPressure::hitPos(
/*WallTriangle* wallTriangle, */double dt, const Particle* p, point_t &hit_pos,
		const point_t &force) {
	hit_pos = p->r + dt * p->v + dt * dt / 2 * force / m_mass;
}

#endif

