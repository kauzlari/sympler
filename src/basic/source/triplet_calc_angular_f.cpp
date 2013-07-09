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
#include "triplet_calc_angular_f.h"

#include "triplet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

const SymbolRegister<TripletCalcAngularF> triplet_calc_angular_f("TripletCalcAngularF");

TripletCalcAngularF::TripletCalcAngularF(Simulation *simulation): TripletCalculator(simulation)
{
  // does not depend on other symbols
  m_stage = 0;

  m_datatype = DataFormat::POINT;
  
  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void TripletCalcAngularF::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalcAngularF");

  m_properties.setDescription( 
			      "This TripletCalculator computes a hookian angular force: F=k(cos(theta)-cos(theta_eq)).");
  
  DOUBLEPC(K, m_k, -1,"Spring force constant.");
  DOUBLEPC(thetaEq, m_thetaEq, -1,"equilibrium angle (in degrees).");
  

  m_k = 0;
  m_thetaEq = 0;

}


void TripletCalcAngularF::setup()
{
  TripletCalculator::setup();

  m_cosEq = cos(PI-PI*m_thetaEq/180.);	
}

void TripletCalcAngularF::setupAfterParticleCreation()
{
  TripletCalculator::setupAfterParticleCreation();
}

#ifndef _OPENMP
void TripletCalcAngularF::compute(triplet_t* tr) {
#else
  void TripletCalcAngularF::compute(triplet_t* tr/*, size_t thread_no*/ /*FIXME: paralelise!*/) {
#endif
      // FIXME!!!: currently not parallelised
    point_t b1, b2; // vector
    double c11 = 0;
    double c12 = 0;
    double c22 = 0; // scalar product
//     double abs_b1, abs_b2;
    double cos_a; //the cosine of the angle
    double F;
//     double a;
    point_t force_a, force_b, force_c;
    // b is the central atom
    for (int _i = 0; _i < SPACE_DIMS; _i++) {
      b1[_i] = tr->b->r[_i] - tr->a-> r[_i]; // b1=rb-ra
      b2[_i] = tr->c->r[_i] - tr->b-> r[_i]; // b2=rc-rb
      // periodic BCs
      if(m_periodic) {
	if(b1[_i] > 0.5*m_boxSize[_i]) b1[_i] -= m_boxSize[_i]; 
	if(b1[_i] < -0.5*m_boxSize[_i]) b1[_i] += m_boxSize[_i]; 
	if(b2[_i] > 0.5*m_boxSize[_i]) b2[_i] -= m_boxSize[_i]; 
	if(b2[_i] < -0.5*m_boxSize[_i]) b2[_i] += m_boxSize[_i]; 
      }      

      c11 += b1[_i] * b1[_i]; // (rb-ra)^2 (dot-product)
      c12 += b1[_i] * b2[_i]; // (rb-ra).(rc-rb)
      c22 += b2[_i] * b2[_i]; // (rc-rb)^2 (dot-product)
    }
    double invAbsC11c22 = 1/sqrt(c11*c22);
    

    // !!! this is PI- the bond angle !!! 
    cos_a = c12*invAbsC11c22; // = (rb-ra).(rc-rb)/(|rb-ra|*|rc-rb|)
    
    //        MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->b->mySlot << " cos(theta) = " << cos_a);
    //       MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->b->mySlot << " m_thetaEq = " << m_thetaEq);
    //       MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->b->mySlot << " cosEQ = " << cos(PI*m_thetaEq/180));
    //        MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->b->mySlot << " m_cosEq = " << m_cosEq);
    
    
    F = -m_k*(cos_a-m_cosEq);
    
    // fa = F*((rb-ra).(rc-rb)/(rb-ra)*(rb-ra)-(rc-rb))/(|rb-ra|*|rc-rb|)
    force_a = F*(c12/c11*b1-b2)*invAbsC11c22; 
    //        MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->a->mySlot << " force_a " << force_a);
    force_c = F*(b1-c12/c22*b2)*invAbsC11c22;
    //        MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->c->mySlot << " force_c " << force_c);
    force_b = -1*force_a - force_c;
    //        MSG_DEBUG("TripletCalcAngularF::calculate force", "particle " << p->b->mySlot << " force_b " << force_b);
    tr->a->tag.pointByOffset(m_slots[0]) += force_a;
    tr->b->tag.pointByOffset(m_slots[1]) += force_b;
    tr->c->tag.pointByOffset(m_slots[2]) += force_c;
    
  }

    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    bool TripletCalcAngularF::findStage()
    {

      // currently (2010/05/17) this always returns true
      return TripletCalculator::findStage();
    }
    
    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    bool TripletCalcAngularF::findStage_0()
    {
      // currently (2010/05/17) this always returns true
      return TripletCalculator::findStage_0();
    }


#ifdef _OPENMP
// FIXME: This module is not yet parallelised. If you parallelise, the following function could roughly do what is commented out now
void TripletCalcAngularF::mergeCopies(size_t thread_no) {
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
