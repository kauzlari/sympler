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


#ifndef __PARTICLE_CACHE_ARB_RNG_H
#define __PARTICLE_CACHE_ARB_RNG_H

#include "particle_cache_arbitrary.h"
#include "random.h"
#include "phase.h"

/*!
 * User-defined random normaly distributed scalar symbol for a particle multiplied by 
 * previously computed particle properties (symbols).
 */
class ParticleCacheArbRNG : public ParticleCacheArbitrary
{
  protected:
   
   /*!
   * a random number generator
   */
   RandomNumberGenerator m_rng;

   /*!
    * pointer to the \a Phase object
    */
   Phase* m_phasePointer;

   /*!
    * Holds user-specified minimum number of particles for which the average of the random numbers will be shifted to zero (by preserving the original variance). Otherwise the average will only be zero over an infinite number of time steps.
    */
   int m_plimit;

    /*!
     * Initialise the PropertyList.
     */
    void init();
  
    
  public:
  /*!
   * Constructor
   */
    ParticleCacheArbRNG(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ParticleCacheArbRNG();

      /*!
     * Compute the cache for particle \a p
     * @param p The particle to compute values for
       */
    virtual void computeCacheFor(Particle* p) = 0;
        
  /*!
     * Setup this ParticleCache
   */
    virtual void setup();
};

#endif
