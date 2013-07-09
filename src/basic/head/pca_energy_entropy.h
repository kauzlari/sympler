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



#ifndef __PARTICLE_CACHE_ENERGY_ENTROPY_H
#define __PARTICLE_CAHCE_ENERGY_ENTROPY_H 

// #include "function.h"
#include "function_fixed.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "weighting_function.h"
#include "val_calculator_rho.h"

/*!
 * This is used by \a FGeneric to calculate the per particle energy and entropy.
 * Otherwise, it would have been calculated more than once per particle.
 * Currently, we do not allow this Cache to be used in the \a Simulation 's list
 */
class ParticleCacheEnergyEntropy: public ParticleCache
{
 protected:
  /*!
   * The derivative of the local configurational energy
   * with respect to the local density
   */
  FunctionFixed *m_de_dn;

  /*!
   * The expression of \a m_de_dn when defined indirectly by some other module.
  */
  string m_de_dn_expr;
  
  /*!
   * The derivative of the local entropy
   * with respect to the local density
   */
  FunctionFixed *m_Tds_dn;
  
  /*!
   * The expression of \a m_Tds_dn when defined indirectly by some other module.
   */
  string m_Tds_dn_expr;

  // m_symbolName is used  
//  string m_de_dn_symbol;

  /*!
   * The symbol of \a m_Tds_dn when defined and used directly in the input file
   */
  string m_Tds_dn_symbol;

  /*!
   * Tag offset of the local energy derivative
   */
  size_t m_de_dn_o;

  /*!
   * Tag offset of the local entropy derivative
   */
  size_t m_Tds_dn_o;


// for the next, we use ParticleCache::m_offset now 
  //  /*!
//   * Tag offset of the local density
  //   */
//  size_t m_density_o;

  /*!
  * Name of the density to be used
  */
  string m_densitySymbol;
  
  /*!
   * Tag offset of the internal energy
   */
  size_t m_energy_o;

  /*!
   * The \a ColourPair (both colors identical)
   */
  ColourPair *m_cp;

  /*!
   * The weighting function to use for this calculation
   */
  WeightingFunction *m_wf;


  /*!
   * Do not allow the density to rise above this value
   */
  double m_density_cutoff;

  /*!
  * Initialise the property list
  */
  virtual void init();
  
 public:
  /*!
   * Constructor
   * @param color Color this cache is used for
   * @param density_cutoff Do not allow the density to rise above this value
   * @param cp The \a ColourPair
   * @param wf The weighting function to use for this calculation
   * @param de_dn Derivative of the local configurational energy
   * @param Tds_dn Derivative of the local entropy
   */
  ParticleCacheEnergyEntropy
    (size_t color/*, size_t offset, string symbolName*/, double density_cutoff, ColourPair *cp, WeightingFunction *wf,
     FunctionFixed *de_dn, FunctionFixed *Tds_dn);

  /*!
   * Constructor for a \a Node list.
  */
  ParticleCacheEnergyEntropy(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~ParticleCacheEnergyEntropy();

  /*!
   * Compute the cached value
   */
  virtual void computeCacheFor(Particle* p) {
    double n = p->tag.doubleByOffset(/*m_density_o*/m_offset);
    double e = p->tag.doubleByOffset(m_energy_o);

    if (n > m_density_cutoff)
      n = m_density_cutoff;

    p->tag.doubleByOffset(m_de_dn_o) = (*m_de_dn)(n, e);
    p->tag.doubleByOffset(m_Tds_dn_o) = (*m_Tds_dn)(n, e);
  }

  /*!
   * Register the appropriate degrees of freedom and depending calculators
   */
  virtual void registerWithParticle();

  /*!
   * Is this calculator equal to \a c
   * @param c Other calculator
   */
  virtual bool operator==(const ParticleCache &c) const {
    if (typeid(c) == typeid(*this)) {
      ParticleCacheEnergyEntropy *cc = (ParticleCacheEnergyEntropy*) &c;

      /* FIXME: expression can be != expression and mean the same! */

      return
	(m_de_dn->expression() == cc->m_de_dn->expression()) &&
          (m_Tds_dn->expression() == cc->m_Tds_dn->expression()) &&
          m_offset == cc->m_offset && m_stage == cc->m_stage && m_colour == cc->m_colour && m_symbolName == cc->m_symbolName;
    } else
      return false;
  }

  /*!
   * Return the cached derivative of the local energy for particle \a p
   * @param p Particle
   */
  double de_dn(Particle *p) const {
    return p->tag.doubleByOffset(m_de_dn_o);
  }

  /*!
   * Return the cached derivative of the local entropy for particle \a p
   * @param p Particle
   */
  double Tds_dn(Particle *p) const {
    return p->tag.doubleByOffset(m_Tds_dn_o);
  }

      /*!
   * Return the name of the computed symbols to be used in other expressions.
       */
  virtual list<string> mySymbolNames()
  {
    list<string> temp;
    assert(temp.empty());
    temp.push_back(m_symbolName);
    temp.push_back(m_Tds_dn_symbol);
    return temp;
  }


  /*!
  * Setup this Calculator
  */
  virtual void setup();
};

#endif
