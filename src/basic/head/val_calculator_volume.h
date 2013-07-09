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



#ifndef __VAL_CALCULATOR_VOLUME_H
#define __VAL_CALCULATOR_VOLUME_H

#include "general.h"
#include "simulation.h"
#include "manager_cell.h"
#include "val_calculator.h"
#include "weighting_function.h"
#include "pca_volume_self_contribution.h"
#include "val_calculator_rho.h"

#define VAL_V_STR "local_volume"



/*!
 * Calculates the local volume
 */
class ValCalculatorVolume : public NonBondedPairParticleCalculator
{
  protected:
  /*!
   * Cut-off radius for local volume calculation. Taken from
   * the weighting function.
   */
    double m_cutoff;

//  /*!
//     * The \a ColourPair for which to calculate the local volume
    //   */
//    ColourPair *m_cp;

  /*!
     * The weighting function to be used for the local volume calculation
   */
    WeightingFunction *m_wf;

  /*!
     * Name of the weighting function
   */
    string m_wfName;

  /*!
  * The tag offset of the local density
    */
    pair<size_t, size_t> m_density_offset;

        /*!
     * The name of the scalar for the local density used by this calculator.
     * If this Calculator is created by other code and if the
     * string is 'none' The calculator creates a ValCalculatorRho to compute it by
     * interpolation. Then this member is set when calling the constructor.
     * If this Calculator is created by the user, it searches for the density
     * attribute and reports an error if it doesn't find it
         */
    string m_rhoSymbol;


//    /*!
//     * Extra string identifier to distinguish between different types of
//     * densities (corrected, uncorrected, etc.)
    //   */
//    string m_extra_string;

  /*!
     * Initialise the property list
   */
    virtual void init();


  public:
  /*!
   * Constructor
   * @param wf Weighting function to use for the local density calculation
   */
    ValCalculatorVolume(string densitySymbol, string symbolName, WeightingFunction *wf/*, string extra_string*/);

   /*!
     * Constructor for the \a Node hierarchy
    */
    ValCalculatorVolume(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorVolume() {
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();

  /*!
     * Return the cut-off radius
   */
    virtual double cutoff() {
      return m_cutoff;
    }

#ifdef _OPENMP
    /*!
     * Merge the copies of all threads together
     */
    virtual void mergeCopies(ColourPair* cp, int thread_no);
#endif

  /*!
     * Compute the local volume constribution of pair \a pD
     * @param pD \a Pairdist whose local volume contribution is calculated
   */
#ifndef _OPENMP
    virtual void compute(Pairdist* pD)
#else
    virtual void compute(Pairdist* pD, int thread_no)
#endif
    {
      if(pD->abs() < m_cutoff) {
        Particle* first;
        Particle* second;
        double Wi, Wj;

        first = pD->firstPart();
        second = pD->secondPart();
//         MSG_DEBUG("ValCalculatorVolume::compute","BEFORE: first = " << first << ", second = " << second);
//         MSG_DEBUG("ValCalculatorVolume::compute","BEFORE: first->c = " << first->c << ", second->c = " << second->c);
//         MSG_DEBUG("ValCalculatorVolume::compute","BEFORE: first->rho = " << first->tag.doubleByOffset(m_density_offset.first) << ", second->rho = " << second->tag.doubleByOffset(m_density_offset.second));
//         MSG_DEBUG("ValCalculatorVolume::compute","BEFORE: first->V = " << first->tag.doubleByOffset(m_slots.first) << ", second->V = " << second->tag.doubleByOffset(m_slots.second));
//         MSG_DEBUG("ValCalculatorVolume::compute","m_density_offset = " << m_density_offset.first << ", " << m_density_offset.second);
//         MSG_DEBUG("ValCalculatorVolume::compute","m_slots = " << m_slots.first << ", " << m_slots.second);

        Wi = m_wf->interpolate(pD, first->r);
        assert(first->tag.doubleByOffset(m_density_offset.first) != 0);
        Wi /= first->tag.doubleByOffset(m_density_offset.first);


        Wj = m_wf->interpolate(pD, second->r);
        assert(second->tag.doubleByOffset(m_density_offset.second) != 0);
        Wj /= second->tag.doubleByOffset(m_density_offset.second);

//         assert(first->c == m_cp->firstColour());
//         assert(second->c == m_cp->secondColour());

//       if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:Wi=" << Wi << ", Wj=" << Wj );
//       if(first->c == 1 && second->c == 0) MSG_DEBUG("ValCalculatorRho::compute", "10:Wi=" << Wi << ", Wj=" << Wj );

/*			assert(m_colour2Slot.find(first->c) != m_colour2Slot.end());
        assert(m_colour2Slot.find(second->c) != m_colour2Slot.end());

        first->tag.doubleByOffset(m_colour2Slot[first->c]) += temp;
        second->tag.doubleByOffset(m_colour2Slot[second->c]) += temp;*/

        if(pD->actsOnFirst())
        {
//           if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorVolume::compute", "01:0=TRUE");
#ifdef ENABLE_PTHREADS
        first->lock();
#endif

#ifndef _OPENMP
        first->tag.doubleByOffset(m_slots.first) += Wj;
#else
        (*first->tag.vectorDoubleByOffset(m_copy_slots[thread_no].first))[m_vector_slots.first] += Wj;
#endif

#ifdef ENABLE_PTHREADS
        first->unlock();
#endif
        }

        if(pD->actsOnSecond())
        {
//           if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorVolume::compute", "01:1=TRUE");
#ifdef ENABLE_PTHREADS
        second->lock();
#endif
//       if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:secondbefore=" << second->tag.doubleByOffset(m_slots.second));
#ifndef _OPENMP
        second->tag.doubleByOffset(m_slots.second) += Wi;
#else
        (*second->tag.vectorDoubleByOffset(m_copy_slots[thread_no].second))[m_vector_slots.second] += Wi;
#endif
//      if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "01:secondafter=" << second->tag.doubleByOffset(m_slots.second));

//      if(first->c == 0 && second->c == 1) MSG_DEBUG("ValCalculatorRho::compute", "2ndslot=" << m_slots.second);

#ifdef ENABLE_PTHREADS
        second->unlock();
#endif
        }
/*        MSG_DEBUG("ValCalculatorVolume::compute","AFTER: first->rho = " << first->tag.doubleByOffset(m_density_offset.first) << ", second->rho = " << second->tag.doubleByOffset(m_density_offset.second));
        MSG_DEBUG("ValCalculatorVolume::compute","AFTER: first->V = " << first->tag.doubleByOffset(m_slots.first) << ", second->V = " << second->tag.doubleByOffset(m_slots.second));*/
      }
    }

  /*!
     * Returns a string identifier containig "ValCalculatorVolume_" and the name
     * of the weighting function as defined in the input file.
   */
    virtual string myName() {
      return m_symbolName /*"ValCalculatorVolume_" + m_extra_string + "_" + m_wf->name()*/;
    }

  /*!
     * Register all degrees of freedom (i.e., the local density)
   */
    virtual void setSlots(ColourPair* cPair, pair<size_t, size_t> &theSlots, bool oneProp) {
      // see CONVENTION5 for rule about persistencies
//       pair<bool, bool> persist(false, false);

//       m_cp = cp;

      cPair->setCutoff(m_cutoff);
      cPair->setNeedPairs(true);
      MSG_DEBUG("ValCalculatorVolume::setSlots", "creating ValCalculatorRho now");
        // first we register the value in the CPs != m_cp, if m_oneProp = true
      if(oneProp)
        FOR_EACH_COLOUR_PAIR
	  (
	   cPair->manager(),
	   // if this is m_cp then do nothing (will be done afterwards)
	   if(cPair->firstColour() != cp->firstColour() || cPair->secondColour() != cp->secondColour())
	     {
	       cp->registerCalc(m_density_offset, new ValCalculatorRho(m_rhoSymbol, m_wf/*, m_extra_string*/), true);
	     }
	   );
            cPair->registerCalc(m_density_offset, new ValCalculatorRho(m_rhoSymbol, m_wf/*, m_extra_string*/), oneProp);
      MSG_DEBUG("ValCalculatorVolume::setSlots", "m_density_offset = (" << m_density_offset.first << ", " << m_density_offset.second << ")");

#define M_SIMULATION ((Simulation*) cPair->manager()->phase()->parent())
#define M_CONTROLLER M_SIMULATION->controller()

// see CONVENTION5 for rule about persistencies
#if 0
// new rules
// FIXME: not checked how meaningful this is for the case that both species don't have an IntegratorPosition
// by the way, how meaningful is this case itself ?!?
      if (!M_CONTROLLER->findIntegrator("IntegratorPosition", cPair->firstSpecies()) &&
        !M_CONTROLLER->findIntegrator("IntegratorPosition", cPair->secondSpecies()))
{
        persist.first = true;
        persist.second = true;
}
#endif

#undef M_CONTROLLER

    MSG_DEBUG
      ("ValCalculatorVolume::setSlots",
       "Registering degree of freedom, cut-off = " << m_cutoff);

    // because of CONVENTION_colour_pair_1, the attributes are only added, if the
    // ValCalculator is newly created

    if(Particle::s_tag_format[cPair->firstColour()].attrExists(m_symbolName))
    {
      if(Particle::s_tag_format[cPair->firstColour()].attrByName(m_symbolName).datatype == DataFormat::DOUBLE)
        m_slots.first =
            Particle::s_tag_format[cPair->firstColour()]./*indexOf*/offsetByName(m_symbolName);
      else
        throw gError("ValCalculatorVolume::setSlots", "Symbol + " + m_symbolName + " already existing as a non-scalar.");
    }
    else
      // see CONVENTION5 for rule about persistencies
      m_slots.first =
        Particle::s_tag_format[cPair->firstColour()].addAttribute
        (/*string(VAL_V_STR) + "_" + m_extra_string + colourString + "_" + m_wf->name()*/m_symbolName,
         DataFormat::DOUBLE, /*persist.first*/false, m_symbolName/*"V[" + m_wf->name() + "]"*/).offset;

    Particle::registerCache
        (new ParticleCacheVolumeSelfContribution
        (cPair->firstColour(), m_slots.first, m_symbolName, m_density_offset.first, m_wf));


    if(cPair->firstColour() == cPair->secondColour()) {
      m_slots.second = m_slots.first;
    } else {
      if(Particle::s_tag_format[cPair->secondColour()].attrExists(m_symbolName))
      {
        if(Particle::s_tag_format[cPair->secondColour()].attrByName(m_symbolName).datatype == DataFormat::DOUBLE)
          m_slots.second =
              Particle::s_tag_format[cPair->secondColour()]./*indexOf*/offsetByName(m_symbolName);
        else
          throw gError("ValCalculatorVolume::setSlots", "Symbol + " + m_symbolName + " already existing as a non-scalar.");
      }
      else
        // see CONVENTION5 for rule about persistencies
        m_slots.second =
            Particle::s_tag_format[cPair->secondColour()].addAttribute
            (/*string(VAL_V_STR) + "_" + m_extra_string + colourString + "_" + m_wf->name()*/m_symbolName,
             DataFormat::DOUBLE, /*persist.second*/false, m_symbolName/*"V[" + m_wf->name() + "]"*/).offset;

      Particle::registerCache
          (new ParticleCacheVolumeSelfContribution
          (cPair->secondColour(), m_slots.second, m_symbolName, m_density_offset.second, m_wf));
    }

    theSlots = m_slots;
#undef M_SIMULATION
    }
};

#endif
