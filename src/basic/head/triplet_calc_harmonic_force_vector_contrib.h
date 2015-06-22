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


#ifndef __TRIPLET_CALC_HARMONIC_FORCE_VECTOR_CONTRIB_H
#define __TRIPLET_CALC_HARMONIC_FORCE_VECTOR_CONTRIB_H

#include "triplet_calculator.h"
#include "triplet.h"

class Pairdist;
class ColourPair;

/*!
 * \a TripletCalculator for cached properties computed during a loop over bonded triplets.
 * This \a Symbol computes the vector part of the force on triplet particles derived from the potential energy function 0.5*psi^2, where psi=pi-alpha, and cos(psi)=cos(pi-alpha)=-cos(alpha)=rij.rjk/(|rij||rjk|); i, j, k are, respectively, the \"left\", \"central\", and \"right\" particle of the triplet, and the vectors rij=ri-rj and rjk=rj-rk. The vector-contributions stored in the particles i, j, k are Fi=rjk/|rij|/|rjk|-cos(psi)*rij/|rij|^2, Fk=cos(psi)*rjk/|rjk|^2-rij/|rij|/|rjk|, and Fj=-Fi-Fk. Currently, this module is limited to the implicitly assumed equilibrium angle psiEq=0 or alphaEq=pi. By default, attribute 'symbol' denotes the variable for storing the value for each particle in the triplet. With attributes 'symbol1st' and 'symbol3rd' this can be modified for the two outer (1st, 3rd) particles. 
 */
class TripletCalcHarmonicForceVectorContrib : public TripletCalculator
{
protected:
  
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
  virtual TripletCalcHarmonicForceVectorContrib* copyMySelf()
  {
    return new TripletCalcHarmonicForceVectorContrib(*this);
  }

  /*!
   * Called by setup() 
   */
  virtual void setupTags();


public:
  /*!
 * Constructor
   */
  TripletCalcHarmonicForceVectorContrib(string symbol);
  
  /*!
   * Constructor
   */
  TripletCalcHarmonicForceVectorContrib(/*Node*/Simulation* parent);
   
  /*!
   * Destructor
   */
  virtual ~TripletCalcHarmonicForceVectorContrib() {
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
/*     double cos_psi; //the cosine of the angle */
    for (int _i = 0; _i < SPACE_DIMS; _i++)
      {
	// -rij=-(ri-rj)
	b1[_i] = tr->b->r[_i] - tr->a-> r[_i];
	// -rjk=-(rj-rk)
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
    
/*     cos_psi = c12*invAbsC11c22; */

    point_t fi /* = cos_psi*b1/c11 - b2*invAbsC11c22 */ = invAbsC11c22*(c12*b1/c11 - b2);
    point_t fk /* = b1*invAbsC11c22 - cos_psi*b2/c22 */ = invAbsC11c22*(b1 - c12*b2/c22);

    // DEBUGGING
/*     if(tr->b->r.y > 19.0) { */
/*       MSG_DEBUG("TripletCalcAllPart::compute", "Triplet debug info:\n" << "p1.r=" << tr->a->r << "\n" << "p2.r=" << tr->b->r << "\n" << "p3.r=" << tr->c->r << "\n"  << "c11=" << c11 << "\n" << "c12=" << c12 << "\n" << "c22=" << c22 << "\n" << "p1.tag before =" <<  tr->a->tag.pointByOffset(m_slots[0]) << "\n" << "p2.tag before =" <<  tr->b->tag.pointByOffset(m_slots[1]) << "\n" << "p3.tag before =" <<  tr->c->tag.pointByOffset(m_slots[2]) << "\n"); */
/*     } */

    tr->a->tag.pointByOffset(m_slots[0]) += fi;
    tr->b->tag.pointByOffset(m_slots[1]) -= (fi+fk);
    tr->c->tag.pointByOffset(m_slots[2]) += fk;

    // DEBUGGING
/*     if(tr->b->r.y > 19.0) { */
/*       MSG_DEBUG("TripletCalcAllPart::compute", "p1.tag after =" <<  tr->a->tag.pointByOffset(m_slots[0]) << "\n" << "p2.tag after =" <<  tr->b->tag.pointByOffset(m_slots[1]) << "\n" << "p3.tag after =" <<  tr->c->tag.pointByOffset(m_slots[2]) << "\n"); */
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
