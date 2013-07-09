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


#ifndef __VAL_CALCULATOR_DIRICHLET_BC_VELS_H
#define __VAL_CALCULATOR_DIRICHLET_BC_VELS_H

#include "val_calculator_BC.h"

#include "particle.h"
#include "pairdist.h"

// #include "wall.h"

/* class Wall; */

//---- newclass -----------------------------------------------------------

/*!
 * Saves the pair-specific value of the velocity of the boundary particle used
 * for applying a Dirichlet boundary condition (BC) in each pair of particles.
 */
class ValCalculatorDirichletBCVels : public ValCalculatorBC
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
      return new ValCalculatorDirichletBCVels(*this); 
    }
  
  public:
  /*!
   * Constructor
   * @param symbol The symbol name
   */
    ValCalculatorDirichletBCVels(string symbol);

  /*!
     * Constructor for Node hierarchy
   */
    ValCalculatorDirichletBCVels(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorDirichletBCVels() {
    }
    
    /*!
     * Register this calculator and save the offset for the data stored
     * in a \a Pairdist in the argument \a slot
     */
/*     void setSlot(ColourPair* cp, size_t& slot, bool oneProp); */
    
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
  
  //   // it does not make sense to check for a cutoff because we should make sure that also the interaction with the longest range may safely use the boundary values computed here
  // /*  if (m_cutoff > pD->abs()) 
  
  Particle* p1st = pD->firstPart();
  Particle* p2nd = pD->secondPart();
  // //          point_t hitpos;
  //     innerDist = HUGE_VAL;
  
  if(m_wallIsSecond)
    {
      //           // if all worked fine, we may now compute the value of the outer particle
      //           // the check for innerDist != HUGE_VAL is done below
      pD->tag.pointByOffset(m_slot) =
	(outerDist+innerDist)*(p2nd->v)/innerDist-(outerDist/innerDist)*p1st->v;
      
/*                  MSG_DEBUG("ValCalculatorDirichletBCVels::compute", "WALLISSECOND: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + ", c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << "r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist); */
      
    }
  else // so !m_wallIsSecond
    {
      
      // if all worked fine, we may now compute the value of the outer particle
      // the check for innerDist != HUGE_VAL is done below
      // the first term is for Dirichlet != 0. We assume the value was assigned 
      // to the wall particles.
      pD->tag.pointByOffset(m_slot) =
	(outerDist+innerDist)*(p1st->v)/innerDist-(outerDist/innerDist)*p2nd->v;

/*       MSG_DEBUG("ValCalculatorDirichletBCVels::compute", "WALLISFIRST: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + ", c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist); */
      
    } // of else of if(m_wallIsSecond)
  
  if(innerDist == HUGE_VAL)
    throw gError("ValCalculatorDirichletBCVels::compute", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + "c1=" + ObjToString(pD->firstPart()->c) + ", c2=" + ObjToString(pD->secondPart()->c) + ", r1=" + ObjToString(pD->firstPart()->r) + ", r2=" + ObjToString(pD->secondPart()->r));
  
  //     //   } // this is the closing of the commented out if(m_cutoff > ...)                                                                          
}

  /*!
     * Returns "dirichletBCVels"
   */
    virtual string myName() {
      return "dirichletBCVels";
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();

    /*!
     * For each wall-particle, find the walls in range
     */
/*     virtual void setupAfterParticleCreation(); */
};

#endif
