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


#ifndef __PCA_RANDOM_H
#define __PCA_RANDOM_H

#include "particle.h"
#include "symbol.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// #define PCA_MAX_STAGE 2

/*!
 * Functions to cache properties for the particles. Similar to the ValCalculators for
 * the pairs.
 */

class PCaRandom: public ParticleCache 
{
 protected:
  /*!
   * Random number distribution.
   */
  string m_distribution;
  
  /*!
   * Standard deviation for gauss distribution.
   */
  double m_sigma;
  
  /*!
   * Random number type.
   */
  const gsl_rng_type *m_rng_type;

  /*!
   * Random number generator.
   */
  gsl_rng *m_rgen;
   
  /*!
   * Seed for the random number generator. If 0, system time is used.
   */
  size_t m_seed;

  /*!
   * Type of the symbol. Can be scalar, vector or tensor.
   */
  string m_type;
  
 
  void init();

  void setup();
  
 public:
  /*!
   * Constructor for Node hierarchy
   */
  PCaRandom(Simulation* parent);
    
  /*!
   * Destructor
   */
  virtual ~PCaRandom();

  /*!
   * Compute the cache for particle \a p
   * @param p The particle to compute values for
   */
  virtual void computeCacheFor(Particle* p);

  /*!
   * Register the additional degrees of freedom with the \a Particle s \a DataFormat
   */
  virtual void registerWithParticle();

  /*!
   * Are those two caches identical?
   * @param c Is this the same cache?
   */

  virtual bool operator==(const ParticleCache &c) const {
    if (typeid(c) == typeid(*this)) {
      PCaRandom *cc = (PCaRandom*) &c;
      return
      m_offset == cc->m_offset && 
      m_symbolName == cc->m_symbolName && 
      m_stage == cc->m_stage && 
      m_colour == cc->m_colour;
    } else
      return false;
  }

  /*!
   * Return the color this cache is for
   */
  size_t colour() const {
    return m_colour;
  }
  
  /*!
   * returns the offset
   */
  size_t offset() const
  {
    return m_offset; 
  }
};
#endif
