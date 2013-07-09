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


#ifndef __BONDED_PAIR_PARTICLE_TENSOR_H
#define __BONDED_PAIR_PARTICLE_TENSOR_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
/* #include "val_calculator_arbitrary.h" */
#include "bonded_pair_particle_arbitrary.h"
#include "colour_pair.h"
#include "function_pair.h"


/*!
 * Function to compute completely user-defined tensor properties for the 
 * particles, which need summation over bonded pairs.
 * FIXME!: This code is unfinished and not yet working as is! 
 */

class BondedPairParticleTensor : public BondedPairParticleArbitrary
/* class BondedPairParticleTensor : public ValCalculatorArbitrary */
{
  protected:

    /*!
     * Initialise the property list
     */
    virtual void init();
    
    /*!
     * Helper function for polymorphic copying
     */
    virtual ValCalculator* copyMySelf()
    {
      return new BondedPairParticleTensor(*this);
    }

    /*!
     * copies the members of this class to \a vc
     */
    virtual void copyMembersTo(ValCalculator* vc);

  
  public:
   
    /*!
   * Constructor for the \a Node hierarchy
     */
    BondedPairParticleTensor(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~BondedPairParticleTensor();

#ifdef _OPENMP
    /*!
     * Merge the copies of all threads together
     */
    virtual void mergeCopies(ColourPair* cp, int thread_no);
#endif

    /*!
     * Compute the user defined expression for pair \a pD
     * @param pD \a Pairdist whose contribution we calculate
     */
#ifndef _OPENMP
        virtual void compute(Pairdist* pD) 
#else
        virtual void compute(Pairdist* pD, int thread_no)
#endif
        {
           
	  point_t temp;
          
	  // compute the pair-expression
	  m_function(&temp, pD);
          
	  tensor_t tempFirst;
	  tensor_t tempSecond;
	  // compute the particle-expressions
	  m_1stparticleFactor(&tempFirst, pD);
	  m_2ndparticleFactor(&tempSecond, pD);
          
	  Particle* first = pD->firstPart(); 
	  Particle* second = pD->secondPart();

	  size_t slot1 = first->mySlot;
	  size_t slot2 = second->mySlot;

	  //             MSG_DEBUG("BondedPairParticleTensor::compute", "temp = " << temp);            
	  

	  if(pD->actsOnFirst())
            {
	      /*This is defined as component-wise multiplication*/
                tempFirst *= temp;

	  if(m_toFile) {
	    o.open(m_fileName.c_str(), ios::out | ios::binary);
	    o.write((char*) &(double(slot1)), sizeof(double));
	    o.write((char*) &(double(slot2)), sizeof(double));
	    o.write((char*) &(tempFirst.tensor), 9*sizeof(double));
	    o.close();
	  }

	      
#ifndef _OPENMP

               first->tag.pointByOffset(m_slots.first) += tempFirst;

/* 	      MSG_DEBUG("BondedPairParticleTensor::compute", "tag-data AFTER = " << first->tag.pointByOffset(m_slots.first));             */


#else
              first->tag.pointByOffset(m_slots.first) += tempFirst;

	      // FIXME: parallelise!

	      // this is how the parallel version could look like
/*               for (size_t a = 0; a < SPACE_DIMS; ++a) { */
/*                 (*first->tag.vectorDoubleByOffset(m_copy_slots[thread_no].first))[m_vector_slots.first + a] += tempFirst[a]; */
/*               } */
#endif
	      
	      //     MSG_DEBUG("BondedPairParticleTensor::compute", "AFTER: first->point = " << first->tag.pointByOffset(m_slots.first));            

            }
	  
	  if(pD->actsOnSecond())
            {
	      /*This is defined as component-wise multiplication*/
                tempSecond *= temp;

	  if(m_toFile) {
	    o.open(m_fileName.c_str(), ios::out | ios::binary);
	    o.write((char*) &(double(slot2)), sizeof(double));
	    o.write((char*) &(double(slot1)), sizeof(double));
	    o.write((char*) &(tempSecond.tensor), 9*sizeof(double));
	    o.close();

	      
#ifndef _OPENMP
              second->tag.tensorByOffset(m_slots.second) += m_symmetry*(tempSecond);
#else
	      // FIXME: parallelise
              second->tag.tensorByOffset(m_slots.second) += m_symmetry*(tempSecond);

#endif
	      
            }
	  
        }

	/*!
         * Setup this Calculator
	 */
        virtual void setup();

};

#endif
