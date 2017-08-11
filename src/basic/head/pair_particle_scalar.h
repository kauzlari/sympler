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


#ifndef __PAIR_PARTICLE_SCALAR_H
#define __PAIR_PARTICLE_SCALAR_H

#include "general.h"
#include "val_calculator_arbitrary.h"
#include "colour_pair.h"


/*!
 * Functions to compute completely user-defined scalar properties for the
 * particles, which need pair summation
 */
class PairParticleScalar : public ValCalculatorArbitrary
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
      return new PairParticleScalar(*this);
    }

  public:

    /*!
   * Constructor for the \a Node hierarchy
     */
    PairParticleScalar(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~PairParticleScalar();

#ifdef _OPENMP
    /*!
     * Merge the copies from all threads together
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
      if(pD->abs() < m_cutoff)
      {
        double temp;

        // compute the expression
        m_function(&temp, pD);

        double tempFirst;
        double tempSecond;
            // compute the particle-expresions
        m_1stparticleFactor(&tempFirst, pD);
        m_2ndparticleFactor(&tempSecond, pD);

        Particle* first = pD->firstPart();
        Particle* second = pD->secondPart();

        if(pD->actsOnFirst())
        {

#ifndef _OPENMP
            first->tag.doubleByOffset(m_slots.first) += temp*tempFirst;
#else

            (*first->tag.vectorDoubleByOffset(m_copy_slots[thread_no].first))[m_vector_slots.first] += temp*tempFirst;

#endif

        }

        if(pD->actsOnSecond())
        {

#ifndef _OPENMP
              second->tag.doubleByOffset(m_slots.second) += m_symmetry*temp*tempSecond;
#else

              (*second->tag.vectorDoubleByOffset(m_copy_slots[thread_no].second))[m_vector_slots.second] += m_symmetry*temp*tempSecond;

#endif

        }
      }
    }

  /*!
     * Returns the symbol name as defined in the input file.
   */
    virtual string myName() {
      return m_symbolName;
    }

  /*!
     * Register all degrees of freedom (i.e., the local density)
   */
    virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp)
    {
      throw gError("PairParticleScalar::setSlots", "should not have been called! Contact the programmer.");
    }

  /*!
     * Setup this Calculator
   */
    virtual void setup();

};

#endif
