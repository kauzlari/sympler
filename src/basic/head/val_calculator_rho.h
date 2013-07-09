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



#ifndef __VAL_CALCULATOR_RHO_H
#define __VAL_CALCULATOR_RHO_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
#include "val_calculator.h"
#include "weighting_function.h"
#include "pca_density_self_contribution.h"
#include "colour_pair.h"

/* #include "valgrind/memcheck.h" */

#define VAL_RHO_STR "local_density"



/*!
 * Calculates the local density via pair summation
 */
class ValCalculatorRho : public NonBondedPairParticleCalculator
/* class ValCalculatorRho : public ValCalculatorPart */
{
 protected:
  /*!
   * Cut-off radius for local density calculation. Taken from
   * the weighting function.
   */
  double m_cutoff;

  /*!
   * The \a ColourPair for which to calculate the local density
   */
  ColourPair *m_cp;

  /*!
   * The weighting function to be used for the local density calculation
   */
  WeightingFunction *m_wf;

  /*!
   * Name of the weighting function
   */
  string m_wfName;

  /*!
   * Should the self contribution of the particles be included into the local
   * density calculation?
  */
  bool m_self;

  /*!
   * Initialise the property list
   */
  virtual void init();


 public:
  /*!
  * Constructor
  * @param wf Weighting function to use for the local density calculation
   */
   ValCalculatorRho(string symbolName, WeightingFunction *wf/*, string extraString = ""*/);/*
    : ValCalculatorPart(), m_cutoff(wf->cutoff()), m_wf(wf) {
    MSG_DEBUG("ValCalculatorRho::ValCalculatorRho", "CONSTRUCTOR");
  }*/

   /*!
   * Constructor for the \a Node hierarchy
   */
   ValCalculatorRho(/*Node*/Simulation* parent);

  /*!
   * Destructor
   */
  virtual ~ValCalculatorRho() {
  }

  /*!
   * Setup this Calculator
   */
  virtual void setup();

#ifdef _OPENMP
  /*!
   * Merge the copies from all threads together
   */
  virtual void mergeCopies(ColourPair* cp, int thread_no);
#endif

  /*!
   * Return the cut-off radius
   */
  virtual double cutoff() {
    return m_cutoff;
  }

  /*!
   * Compute the local density constribution of pair \a pD
   * @param pD \a Pairdist whose local density contribution is calculated
   */
#ifndef _OPENMP
  virtual void compute(Pairdist* pD)
#else
  virtual void compute(Pairdist* pD, int thread_no)
#endif
  {
//     MSG_DEBUG("ValCalculatorRho::compute","called:" << pD->firstPart()->c << ", " << pD->secondPart()->c);

    if(pD->abs() < m_cutoff) {
      Particle* first;
      Particle* second;
      double Wi, Wj;

      first = pD->firstPart();
      second = pD->secondPart();

      Wi = m_wf->interpolate(pD, first->r);
      Wj = m_wf->interpolate(pD, second->r);

//       MSG_DEBUG("ValCalculatorRho::compute","Wi = " << Wi << ", Wj = " << Wj << ", slots = " << m_slots.first << "," << m_slots.second);


//       assert(first->c == m_cp->firstColour());
//       assert(second->c == m_cp->secondColour());

//       if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:Wi=" << Wi << ", Wj=" << Wj );
//       if(first->c == 1 && second->c == 0) MSG_DEBUG("ValCalculatorRho::compute", "10:Wi=" << Wi << ", Wj=" << Wj );

/*			assert(m_colour2Slot.find(first->c) != m_colour2Slot.end());
			assert(m_colour2Slot.find(second->c) != m_colour2Slot.end());

			first->tag.doubleByOffset(m_colour2Slot[first->c]) += temp;
			second->tag.doubleByOffset(m_colour2Slot[second->c]) += temp;*/

      if(pD->actsOnFirst())
      {
/*if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:0=TRUE");*/
#ifdef ENABLE_PTHREADS
        first->lock();
#endif

#ifndef _OPENMP
        first->tag.doubleByOffset(m_slots.first) += Wj;
#else

/* 	MSG_DEBUG("ValCalculatorRho::compute", "m_vector_slots.first = " << m_vector_slots.first ); */

        (*first->tag.vectorDoubleByOffset(m_copy_slots[thread_no].first))[m_vector_slots.first] += Wj;
#endif

// MSG_DEBUG("ValCalculatorRho::compute", "firstafter=" << first->tag.doubleByOffset(m_slots.second));
#ifdef ENABLE_PTHREADS
        first->unlock();
#endif
      }

      if(pD->actsOnSecond())
      {
/*if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:1=TRUE");*/
#ifdef ENABLE_PTHREADS
        second->lock();
#endif
//       if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:secondbefore=" << second->tag.doubleByOffset(m_slots.second));
#ifndef _OPENMP
        second->tag.doubleByOffset(m_slots.second) += Wi;
#else
        (*second->tag.vectorDoubleByOffset(m_copy_slots[thread_no].second))[m_vector_slots.second] += Wi;
#endif
// MSG_DEBUG("ValCalculatorRho::compute", "secondafter=" << second->tag.doubleByOffset(m_slots.second));

//      if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "2ndslot=" << m_slots.second);



#ifdef ENABLE_PTHREADS
        second->unlock();
#endif
      }
    }
  }

  /*!
   * Returns the symbol name as defined in the input file.
   */
  virtual string myName() {
    return m_symbolName /*"ValCalculatorRho_" + m_extra_string + "_" + m_wf->name()*/;
  }

  /*!
   * Register all degrees of freedom (i.e., the local density)
   */
  virtual void setSlots(ColourPair* cp, pair<size_t, size_t> &theSlots, bool oneProp) {
    // see CONVENTION5 for rule about persistencies
//     pair<bool, bool> persist(false, false);

//     m_cp = cp;

    /*m_*/cp->setCutoff(m_cutoff);
    /*m_*/cp->setNeedPairs(true);

#define M_SIMULATION ((Simulation*) cp->manager()->phase()->parent())

/*    MSG_DEBUG
      ("ValCalculatorRho::setSlots",
       "Registering degree of freedom, cut-off = " << m_cutoff);*/

string colourString;
    // are the properties ColourPair-UNSPECIFIC?
    if(/*M_SIMULATION->*/oneProp/*()*/)
    {
      colourString = "";
    }
    // the properties are ColourPair-SPECIFIC
    else
    {
      colourString = "_" + cp->toString();
    }

    // because of CONVENTION_colour_pair_1, the attributes are only added, if the
    // ValCalculator is newly created

    // see CONVENTION5 for rule about persistencies
    m_slots.first =
      Particle::s_tag_format[cp->firstColour()].addAttribute
        (m_symbolName, DataFormat::DOUBLE, false, m_symbolName).offset;

    Particle::registerCache
      (new ParticleCacheDensitySelfContribution
       (/*m_*/cp->firstColour(), m_slots.first, m_wf, m_symbolName));


    if(cp->firstColour() == cp->secondColour()) {
      m_slots.second = m_slots.first;
    } else {
      // see CONVENTION5 for rule about persistencies
      m_slots.second =
        Particle::s_tag_format[cp->secondColour()].addAttribute
          (m_symbolName, DataFormat::DOUBLE, false, m_symbolName).offset;

    Particle::registerCache
      (new ParticleCacheDensitySelfContribution
       (/*m_*/cp->secondColour(), m_slots.second, m_wf, m_symbolName));
    }

    theSlots = m_slots;
#undef M_SIMULATION
  }
};

#endif
