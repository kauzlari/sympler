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
#include "triplet_calc_central_part_scalar.h"

#include "triplet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

const SymbolRegister<TripletCalcCentralPartScalar> triplet_calc_central_part_scalar("TripletCalcCentralPartScalar");

TripletCalcCentralPartScalar::TripletCalcCentralPartScalar(Simulation *simulation): TripletCalculator(simulation)
{
  // does not depend on other symbols
  m_stage = 0;

  m_datatype = DataFormat::DOUBLE;
  
  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void TripletCalcCentralPartScalar::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalcCentralPartScalar");

  m_properties.setDescription("TripletCalculator for cached properties computed during a loop over bonded triplets. This TripletCalculator specifically computes a user-defined property depending on the cosine of the triplet angle \"ca\". For each triplet this property only contributes to the CENTRE particle by summation! ");
  
  FUNCTIONFIXEDPC
      (expression, m_expression, "Function depending on the variable \"ca\" which is the cosine of the triplet-angle.");
  m_expression.addVariable("ca");

}


void TripletCalcCentralPartScalar::setup()
{
  TripletCalculator::setup();
}

void TripletCalcCentralPartScalar::setupAfterParticleCreation()
{
  TripletCalculator::setupAfterParticleCreation();
}

#ifndef _OPENMP
void TripletCalcCentralPartScalar::compute(triplet_t* tr) {
#else
  void TripletCalcCentralPartScalar::compute(triplet_t* tr/*, size_t thread_no*/ /*FIXME: paralelise!*/) {
#endif
      // FIXME!!!: currently not parallelised
    point_t b1, b2; // vector
    double c11 = 0;
    double c12 = 0;
    double c22 = 0; // scalar product
    double cos_a; //the cosine of the angle
    for (int _i = 0; _i < SPACE_DIMS; _i++)
      {
	b1[_i] = tr->b->r[_i] - tr->a-> r[_i];
	b2[_i] = tr->c->r[_i] - tr->b-> r[_i];
	// periodic BCs
	if(m_periodic) {
	  if(b1[_i] > 0.5*m_boxSize[_i]) b1[_i] -= m_boxSize[_i]; 
	  if(b1[_i] < -0.5*m_boxSize[_i]) b1[_i] += m_boxSize[_i]; 
	  if(b2[_i] > 0.5*m_boxSize[_i]) b2[_i] -= m_boxSize[_i]; 
	  if(b2[_i] < -0.5*m_boxSize[_i]) b2[_i] += m_boxSize[_i]; 
	}	

	c11 += b1[_i] * b1[_i];
	c12 += b1[_i] * b2[_i];
	c22 += b2[_i] * b2[_i];
      }
    double invAbsC11c22 = 1/sqrt(c11*c22);
    
    // if the desired angle is a, the next WITHOUT THE MINUS gives in fact cos(pi-a)=-cos(a) 
    cos_a = -c12*invAbsC11c22;
    
    //        MSG_DEBUG("TripletCalcCentralPartScalar::calculate force", "particle " << p->b->mySlot << " cos(theta) = " << cos_a);
    //       MSG_DEBUG("TripletCalcCentralPartScalar::calculate force", "particle " << p->b->mySlot << " m_thetaEq = " << m_thetaEq);
    //        MSG_DEBUG("TripletCalcCentralPartScalar::calculate force", "particle " << p->b->mySlot << " m_cosEq = " << m_cosEq);
        
    // contributes only to central particle! (because, currently(2010-05-18), more not needed)
    tr->b->tag.doubleByOffset(m_slots[1]) += m_expression(cos_a);
    
  }

    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    bool TripletCalcCentralPartScalar::findStage()
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
    bool TripletCalcCentralPartScalar::findStage_0()
    {
      // currently (2010/05/17) this always returns true
      return TripletCalculator::findStage_0();
    }


#ifdef _OPENMP
// FIXME: This module is not yet parallelised. If you parallelise, the following function could roughly do what is commented out now
void TripletCalcCentralPartScalar::mergeCopies(size_t thread_no) {
//   size_t slot2 = m_slots[1];

//   size_t copySlot2 = m_copy_slots[thread_no][1];

//   size_t vecSlot2 = m_vector_slots[1];

//   FOR_EACH_PARTICLE_C 
//   (M_PHASE, m_secondColour,
//     for (size_t j = 0; j < SPACE_DIMS; ++j) {
//       __iSLFE->tag.doubleByOffset(slot2)[j] += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j];
//       (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2 + j] = 0;
//     }
//   );
}
#endif
