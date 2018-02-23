/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


#ifndef __PARTICLE_CACHE_DENSITY_0OC_H
#define __PARTICLE_CAHCE_DENSITY_0OC_H

#include "particle_cache.h"
#include "weighting_function.h"

class ParticleCacheDensity0Oc: public ParticleCache
{
  protected:
    
  /*!
   * Tag offset of the local volume integral
   */
    size_t m_volume_offset;

  /*!
     * Name to be given to the Symbol for the local volume of a particle
   */
    string m_volumeSymbolName;
        
//    /*!
//     * Tag offset of the local density
    //   */
//    size_t m_density_offset;

    
    // with the new hierarchy, the next shouldn't be needed, because only the 
    // uncorrected density and the local volume need the weighting function    

     /*!
    * Weighting function to use
      */
   WeightingFunction *m_wf;

  /*!
     * The name of the weigting function.
   */
    string m_weighting_function;

    /*!
    * Initialise the property list
    */
    virtual void init();

    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new ParticleCacheDensity0Oc(*this);
    }
   
 public:
    
    /*!
     * Constructor
     */
    ParticleCacheDensity0Oc
      (/*Node*/Simulation* parent);

    /*!
     * Destructor
     */
    virtual ~ParticleCacheDensity0Oc();

  /*!
     * Correct the particle's density by dividing through its volume integral
   */  
    virtual void computeCacheFor(Particle* p) {
      assert(p->tag.doubleByOffset(m_volume_offset) != 0.);
      p->tag.doubleByOffset(/*m_density_offset*/m_offset) /= p->tag.doubleByOffset(m_volume_offset);
    }

  /*!
     * Take steps necessary to register this calculator
   */
    virtual void registerWithParticle();

  /*!
     * Does this calculator equal \a c?
     * @param c Other calculator
   */
    virtual bool operator==(const ParticleCache &c) const {
      if (typeid(c) == typeid(*this)) {
        ParticleCacheDensity0Oc *cc = (ParticleCacheDensity0Oc*) &c;

        return
            m_wf->name() == cc->m_wf->name() && m_volume_offset == cc->m_volume_offset &&
            /*m_density_offset*/m_offset == cc->/*m_density_offset*/m_offset && m_colour == cc->m_colour && m_stage == cc->m_stage && m_symbolName == cc->m_symbolName;
      } else
        return false;
    }

  /*!
     * Return name of the Symbol for the local volume of a particle
   */
    string myVolumeSymbolName()
    {
      return m_volumeSymbolName; 
    }
    
  /*!
     * Setup this ParticleCacheDensity0Oc
   */
    virtual void setup();

};

#endif
