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
#include "quintet_calc_curvature.h"

#include "quintet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

const SymbolRegister<QuintetCalcCurvature> quintet_calc_curvature("QuintetCalcCurvature");

QuintetCalcCurvature::QuintetCalcCurvature(Simulation *simulation): QuintetCalculator(simulation)
{
  // does not depend on other symbols
  m_stage = 0;

  m_datatype = DataFormat::DOUBLE;
  
  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void QuintetCalcCurvature::init()
{
  m_properties.setClassName("QuintetCalcPart");
  m_properties.setName("QuintetCalcCurvature");

  m_properties.setDescription( 
			      "This QuintetCalculator computes the curvature at the certral particle p11 using a surface represented by five particles. The five particle are assumed to form a surface that can be expressed by lagrange polynomials. The calculator is only valid in R3.");
}

void QuintetCalcCurvature::setup()
{
  QuintetCalculator::setup();

}

void QuintetCalcCurvature::setupAfterParticleCreation()
{
  QuintetCalculator::setupAfterParticleCreation();
}

#ifndef _OPENMP
void QuintetCalcCurvature::compute(quintet_t* q) {
#else
  void QuintetCalcCurvature::compute(quintet_t* q/*, size_t thread_no*/ /*FIXME: paralelise!*/) {
#endif
      point_t b00, b20, b02, b22; // position vectors
      point_t nv;
      point_t Ss, St; 	

      double E = 0;
      double F = 0;
      double G = 0;

      double N = 0;
      double M = 0;
      double L = 0;

      //int _1 = -1;
      int _j = 0;
      int _k = 0;

      for (int _i = 0; _i < SPACE_DIMS; _i++)
	{   
	  b00[_i] = q->p00->r[_i] - q->p11->r[_i];
	  b20[_i] = q->p20->r[_i] - q->p11->r[_i];
	  b02[_i] = q->p02->r[_i] - q->p11->r[_i];
	  b22[_i] = q->p22->r[_i] - q->p11->r[_i];

	  // periodic BCs
	  if(m_periodic == true)
	    {
	      if(b00[_i] > 0.5*m_boxSize[_i]) b00[_i] -= m_boxSize[_i]; 
	      if(b00[_i] < -0.5*m_boxSize[_i]) b00[_i] += m_boxSize[_i]; 
	      if(b20[_i] > 0.5*m_boxSize[_i]) b20[_i] -= m_boxSize[_i]; 
	      if(b20[_i] < -0.5*m_boxSize[_i]) b20[_i] += m_boxSize[_i];
 	      if(b02[_i] > 0.5*m_boxSize[_i]) b02[_i] -= m_boxSize[_i]; 
	      if(b02[_i] < -0.5*m_boxSize[_i]) b02[_i] += m_boxSize[_i]; 
	      if(b22[_i] > 0.5*m_boxSize[_i]) b22[_i] -= m_boxSize[_i]; 
	      if(b22[_i] < -0.5*m_boxSize[_i]) b22[_i] += m_boxSize[_i]; 
	  }

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
	}	

	L = N;

	// The mean curvature caculates as H = A/B*1/2 
	double A = E*N + G*L - 2*F*M;
	double B = E*G-F*F; 

	double H = A/B*0.5;


    q->p00->tag.doubleByOffset(m_slots[0]) += 0;
    q->p02->tag.doubleByOffset(m_slots[1]) += 0;
    q->p22->tag.doubleByOffset(m_slots[2]) += 0;
    q->p20->tag.doubleByOffset(m_slots[3]) += 0;
    q->p11->tag.doubleByOffset(m_slots[4]) = H;
  }

    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    bool QuintetCalcCurvature::findStage()
    {

      // currently (2010/05/17) this always returns true
      return QuintetCalculator::findStage();
    }
    
    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    bool QuintetCalcCurvature::findStage_0()
    {
      // currently (2010/05/17) this always returns true
      return QuintetCalculator::findStage_0();
    }


#ifdef _OPENMP
// FIXME: This module is not yet parallelised. If you parallelise, the following function could roughly do what is commented out now
void QuintetCalcCurvature::mergeCopies(size_t thread_no) {
//   size_t slot1 = m_slots[0];
//   size_t slot2 = m_slots[1];
//   size_t slot3 = m_slots[2];

//   size_t copySlot1 = m_copy_slots[thread_no][0];
//   size_t copySlot2 = m_copy_slots[thread_no][1];
//   size_t copySlot3 = m_copy_slots[thread_no][2];
//   size_t vecSlot1 = m_vector_slots[0];
//   size_t vecSlot2 = m_vector_slots[1];
//   size_t vecSlot3 = m_vector_slots[2];

//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, m_firstColour,
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slot1)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1 + j] = 0;
//     }
//   );
//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, m_secondColour,
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slot2)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j] = 0;
//     }
//   );
//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, m_thirdColour,
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.pointByOffset(slot3)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot3))[vecSlot3 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot3))[vecSlot3 + j] = 0;
//     }
//   );
}
#endif
