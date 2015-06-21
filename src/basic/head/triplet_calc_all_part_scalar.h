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


#ifndef __TRIPLET_CALC_ALL_PART_SCALAR_H
#define __TRIPLET_CALC_ALL_PART_SCALAR_H

#include "triplet_calc_part_scalar.h"
#include "triplet.h"

class Pairdist;
class ColourPair;

/*!
 * \a TripletCalculator for cached properties computed during a loop over bonded triplets.
 * This \a Symbol computes (possibly different) runtme-compiled expressions depending on the
 * cosine of triplet angle "cosa" and sums them up in the respective particles of the triplet. 
 */
class TripletCalcAllPartScalar : public TripletCalcPartScalar
{
protected:
  
  /*!
   * Additional user-defined function for the first particle depending on the cosine of the triplet-angle
   */
  FunctionFixed m_expression1st;
  
  /*!
   * Additional user-defined function for the third particle depending on the cosine of the triplet-angle
   */
  FunctionFixed m_expression3rd;

  /*!
   * Additional symbol name for the first particle of the triplet
   */
  string m_symbolName1st;

  /*!
   * Additional symbol name for the third particle of the triplet
   */
  string m_symbolName3rd;
  
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   */
  virtual TripletCalcAllPartScalar* copyMySelf()
  {
    return new TripletCalcAllPartScalar(*this);
  }

  /*!
   * Called by setup() 
   */
  virtual void setupTags();


public:
  /*!
 * Constructor
   */
  TripletCalcAllPartScalar(string symbol);
  
  /*!
   * Constructor
   */
  TripletCalcAllPartScalar(/*Node*/Simulation* parent);
   
  /*!
   * Destructor
   */
  virtual ~TripletCalcAllPartScalar() {
  }

  /*!
   * Compute cached properties for \a triplet_t \a tr
   * @param tr The \a triplet_t for which to compute the cached properties
   */
#ifndef _OPENMP
  virtual void compute(triplet_t* tr) {
#else
    virtual void compute(triplet_t* tr/*, size_t thread_no*/ /*FIXME: paralelise!*/) {
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

    // DEBUGGING
/*     if(tr->b->r.y > 19.0) { */
/*       MSG_DEBUG("TripletCalcAllPart::compute", "Triplet debug info:\n" << "p1.r=" << tr->a->r << "\n" << "p2.r=" << tr->b->r << "\n" << "p3.r=" << tr->c->r << "\n"  << "c11=" << c11 << "\n" << "c12=" << c12 << "\n" << "c22=" << c22 << "\n"  << "cosa=" << cos_a << "\n"  << "alpha=" << acos(cos_a) << "\n" << "p1.tag before =" <<  tr->a->tag.doubleByOffset(m_slots[0]) << "\n" << "expressionL=" << m_expression1st.expression() << "\n" << "expressionC=" << m_expression.expression() << "\n" << "expressionR=" << m_expression3rd.expression() << "\n" << "p2.tag before =" <<  tr->b->tag.doubleByOffset(m_slots[1]) << "\n" << "p3.tag before =" <<  tr->c->tag.doubleByOffset(m_slots[2]) << "\n"); */
/*     } */

    tr->a->tag.doubleByOffset(m_slots[0]) += m_expression1st(cos_a);
    tr->b->tag.doubleByOffset(m_slots[1]) += m_expression(cos_a);
    tr->c->tag.doubleByOffset(m_slots[2]) += m_expression3rd(cos_a);

    // DEBUGGING
/*     if(tr->b->r.y > 19.0) { */
/*       MSG_DEBUG("TripletCalcAllPart::compute", "p1.tag after =" <<  tr->a->tag.doubleByOffset(m_slots[0]) << "\n" << "p2.tag after =" <<  tr->b->tag.doubleByOffset(m_slots[1]) << "\n" << "p3.tag after =" <<  tr->c->tag.doubleByOffset(m_slots[2]) << "\n"); */
/*     } */

    
  }


#ifdef _OPENMP
  /*!
   * Merging the copies among different threads (processors) together
   */
    virtual void mergeCopies(size_t thread_no);
#endif

    /*!
     * setup this calculator
     */
  virtual void setup();

    /*!
     * run those setups requiring that \a Particle s have been created
     */
  virtual void setupAfterParticleCreation();


    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
  bool findStage();

    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
  bool findStage_0();


};


#endif
