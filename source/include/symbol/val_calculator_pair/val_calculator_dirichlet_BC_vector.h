/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018,
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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

#ifndef __VAL_CALCULATOR_DIRICHLET_BC_VECTOR_H
#define __VAL_CALCULATOR_DIRICHLET_BC_VECTOR_H

#include "val_calculator_dirichlet_BC_arbitrary.h"

/*!
 * Saves the pair-specific value of an arbitrary vector of the boundary
 * particle used for applying a Dirichlet boundary condition (BC) in each pair
 * of particles.
 */
class ValCalculatorDirichletBCVector : public ValCalculatorDirichletBCArbitrary
{

	protected:

		/*!
     * Initialise the property list
     */
    virtual void init();
  
    /*!
     * Helper function for returning a copy of itself
     */
    virtual ValCalculatorPair* copyMySelf()
    {
      return new ValCalculatorDirichletBCVector(*this);
    }

  public:

    /*!
     * Constructor
     * @param symbol The symbol name
     */
    ValCalculatorDirichletBCVector(string symbol);

    /*!
     * Constructor for Node hierarchy
     */
    ValCalculatorDirichletBCVector(/*Node*/Simulation* parent);

    /*!
     * Destructor
     */
    virtual ~ValCalculatorDirichletBCVector() {
    }

    /*!
     * Compute the boundary value for the \a Pairdist \a pD
     * @param pD \a Pairdist for which to compute the boundary value
     */ 
#ifndef _OPENMP
    virtual void compute(Pairdist* pD)
#else
      virtual void compute(Pairdist* pD, int thread_no)
#endif
    {
      double innerDist;
      double outerDist;
      
      computeDists(pD, innerDist, outerDist);
      
      Particle* p1st = pD->firstPart();
      Particle* p2nd = pD->secondPart();
      
      if(m_wallIsSecond) {
	  
      	// if all worked fine, we may now compute the value of the outer particle
      	// the check for innerDist != HUGE_VAL is done below
      	pD->tag.pointByOffset(m_slot) =
      			(outerDist+innerDist)
						* (p2nd->tag.pointByOffset(m_varOffset.second)) / innerDist
						- (outerDist / innerDist)
							* p1st->tag.pointByOffset(m_varOffset.first);
	  
      }
      else { // so !m_wallIsSecond
	  
      	// if all worked fine, we may now compute the value of the outer particle
      	// the check for innerDist != HUGE_VAL is done below
      	// the first term is for Dirichlet != 0. We assume the value was assigned
      	// to the wall particles.
      	pD->tag.pointByOffset(m_slot) =
      			(outerDist+innerDist)
						* (p1st->tag.pointByOffset(m_varOffset.first)) / innerDist
						- (outerDist / innerDist)
							* p2nd->tag.pointByOffset(m_varOffset.second);
	  
      } // of else of if(m_wallIsSecond)
      
      if(innerDist == HUGE_VAL)
      	throw gError("ValCalculatorDirichletBCVector::compute", "No wall found "
      			"for pair. Check your geometry and other settings. If this doesn't "
      			"help, contact the programmers. \nDetails: slot1="
      			+ ObjToString(pD->firstPart()->mySlot) + ", slot2="
						+ ObjToString(pD->secondPart()->mySlot) + "c1="
						+ ObjToString(pD->firstPart()->c) + ", c2="
						+ ObjToString(pD->secondPart()->c) + ", r1="
						+ ObjToString(pD->firstPart()->r) + ", r2="
						+ ObjToString(pD->secondPart()->r));
    }

    /*!
     * Returns "dirichletBCVector"
     */
    virtual string myName() {
      return "dirichletBCVector";
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();

};

#endif
