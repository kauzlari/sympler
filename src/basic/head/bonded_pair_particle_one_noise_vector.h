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


#ifndef __BONDED_PAIR_PARTICLE_ONE_NOISE_VECTOR_H
#define __BONDED_PAIR_PARTICLE_ONE_NOISE_VECTOR_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
#include "bonded_pair_particle_calc.h"
#include "colour_pair.h"
#include "function_pair.h"


/*!
 * Class doing summation over bonded pairs and computing a DPD-like momentum conserving vector force with random amplitude. 
 *
 *ASSUMPTION: constant integration time step
 *
 */
class BondedPairParticleOneNoiseVector : public BondedPairParticleCalc
{
 protected:

  RandomNumberGenerator m_rng;

  bool m_randomize;

  /*!
   * The noise amplitude
   */
  double m_noise;
  
  /*!
   * = m_noise/sqrt(dt)
   */
  double m_noise_and_time;


  
  /*!
   * Initialise the property list
   */
  virtual void init();
  
  /*!
   * Helper function for polymorphic copying
   */
  virtual ValCalculator* copyMySelf()
  {
    return new BondedPairParticleOneNoiseVector(*this);
  }
  
  public:
   
    /*!
   * Constructor for the \a Node hierarchy
     */
    BondedPairParticleOneNoiseVector(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~BondedPairParticleOneNoiseVector();

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

	  point_t temp = pD->cartesian()/pD->abs();

          temp*=m_rng.normal(1) * m_noise_and_time;
                    
	  Particle* first = pD->firstPart(); 
	  Particle* second = pD->secondPart();
          
	  if(pD->actsOnFirst())
            {
	      
#ifndef _OPENMP

               first->tag.pointByOffset(m_slots.first) += temp;

#else
	      // FIXME: parallelise!
              first->tag.pointByOffset(m_slots.first) += temp;
#endif
	      
            }
	  
	  if(pD->actsOnSecond())
            {
	      
#ifndef _OPENMP
              second->tag.pointByOffset(m_slots.second) += m_symmetry*temp;
#else
	      // FIXME: parallelise
              second->tag.pointByOffset(m_slots.second) += m_symmetry*temp;
#endif
	      
            }
	  
        }

	/*!
         * Returns the symbol name as defined in the input file.
	 */
        virtual string myName() {
          return m_symbolName;
        }

	/*!
         * Returns the name of the bonded list this \a ValCalculator computes for.
	 */
        virtual string listName() {
          return m_listName;
        }

	/*!
         * Register the computed symbol
	 */
        virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) 
        {
          throw gError("BondedPairParticleOneNoiseVector::setSlots", "should not have been called! Contact the programmer.");
        }

	/*!
         * Setup this Calculator
	 */
        virtual void setup();    

	/*!
	 * copies the members of this class to \a vc
	 */
	virtual void copyMembersTo(ValCalculator* vc);



};

#endif
