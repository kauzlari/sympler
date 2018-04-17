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


#ifndef __PAIR_PARTICLE_TENSOR_H
#define __PAIR_PARTICLE_TENSOR_H

#include "general.h"
// #include "simulation.h"
// #include "manager_cell.h"
#include "val_calculator_arbitrary.h"
// #include "colour_pair.h"


/*!
 * Function to compute completely user-defined tensor properties for the
 * particles, which need pair summation
 */
class PairParticleTensor : public ValCalculatorArbitrary
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
      return new PairParticleTensor(*this);
    }

  public:

    /*!
   * Constructor for the \a Node hierarchy
     */
    PairParticleTensor(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~PairParticleTensor();

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
//           MSG_DEBUG("PairParticleVector::compute", "START, m_cutoff = " << m_cutoff << ", slots = (" << m_slots.first << ", " << m_slots.second << ")");
      if(pD->abs() < m_cutoff)
      {
        tensor_t temp;

            // compute the pair-expression
        m_function(&temp, pD);

        tensor_t tempFirst;
        tensor_t tempSecond;
            // compute the particle-expressions
        m_1stparticleFactor(&tempFirst, pD);
        m_2ndparticleFactor(&tempSecond, pD);

        Particle* first = pD->firstPart();
        Particle* second = pD->secondPart();

        if(pD->actsOnFirst())
        {
          for(size_t i = 0; i < SPACE_DIMS; ++i)
            for(size_t j = 0; j < SPACE_DIMS; ++j)
              tempFirst(i, j) *= temp(i, j);

#ifndef _OPENMP
              first->tag.tensorByOffset(m_slots.first) += tempFirst;
#else
              size_t _i = 0;
              for(size_t a = 0; a < SPACE_DIMS; ++a) {
                for(size_t b = 0; b < SPACE_DIMS; ++b) {
                  (*first->tag.vectorDoubleByOffset(m_copy_slots[thread_no].first))[m_vector_slots.first + _i] += tempFirst(a, b);
                  ++_i;
                }
              }
#endif

        }

        if(pD->actsOnSecond())
        {
          for(size_t i = 0; i < SPACE_DIMS; ++i)
            for(size_t j = 0; j < SPACE_DIMS; ++j)
              tempSecond(i, j) *= temp(i, j);

#ifndef _OPENMP
              second->tag.tensorByOffset(m_slots.second) += m_symmetry*(tempSecond);
#else
              size_t _i = 0;
              for(size_t a = 0; a < SPACE_DIMS; ++a) {
                for(size_t b = 0; b < SPACE_DIMS; ++b) {
                  (*second->tag.vectorDoubleByOffset(m_copy_slots[thread_no].second))[m_vector_slots.second + _i] += m_symmetry*tempSecond(a, b);
                  ++_i;
                }
              }
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
#ifdef _OPENMP
   /*!
    *
    */
    virtual void mergeCopies(ColourPair* cp, int thread_no);
#endif

  /*!
     * Register all degrees of freedom (i.e., the local density)
   */
    virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp)
    {
      throw gError("PairParticleTensor::setSlots", "should not have been called! Contact the programmer.");
    }

  /*!
     * Setup this Calculator
   */
    virtual void setup();

};

#endif
