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


#ifndef __PARTICLE_CACHE_VOLUME_SELF_CONTRIBUTION_H
#define __PARTICLE_CAHCE_VOLUME_SELF_CONTRIBUTION_H 

#include "particle_cache.h"
#include "weighting_function.h"

/*!
 * Adds the particles self contribution to the local volume
 * to the local volume field.
 */
class ParticleCacheVolumeSelfContribution: public ParticleCache
{
  protected:
  
//    /*!
//   * Tag offset of the local volume
    //   */
//    size_t m_offset;

  /*!
     * Tag offset of the local density
   */
    size_t m_density_offset;

    /*!
    * From which symbol to take the local density
    */
    string m_densitySymbol;
    
  /*!
     * Weighting function to use for the local volume calculation
   */
    WeightingFunction *m_wf;
    
    /*!
    * Name of the weighting function
    */
    string m_wfName;
    
    /*!
    * Initialise the property list
    */
    virtual void init();

    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new ParticleCacheVolumeSelfContribution(*this);
    }

    
  public:
  /*!
   * Constructor
   * @param color The particle's color
   * @param offset Tag offset of the local volume
   * @param density_offset Tag offset of the local density
   * @param wf The weighting function to use
   */
    ParticleCacheVolumeSelfContribution
        (size_t color, size_t offset, string symbolName, size_t density_offset, WeightingFunction *wf);
    
    /*!
    * Constructor for node hierarchy
    */
    ParticleCacheVolumeSelfContribution(/*Node*/Simulation* parent);
 
  /*!
     * Destructor
   */
    virtual ~ParticleCacheVolumeSelfContribution();

  /*!
     * Compute the self contribution to the local volume
   */
    virtual void computeCacheFor(Particle* p) {
      p->tag.doubleByOffset(m_offset) += m_wf->interpolate(NULL, p->r) / p->tag.doubleByOffset(m_density_offset);
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
      MSG_DEBUG("ParticleCacheVolumeSelfContribution::==", "called");
      if (typeid(c) == typeid(*this)) {
        ParticleCacheVolumeSelfContribution *cc = (ParticleCacheVolumeSelfContribution*) &c;

        return
            m_wf->name() == cc->m_wf->name() && m_offset == cc->m_offset && m_density_offset == cc->m_density_offset && m_stage == cc->m_stage && m_offset == cc->m_offset && m_colour == cc->m_colour;
      } else
        return false;
    }

    /*!
    * Setup this Calculator
    */
    virtual void setup();
    
};

#endif
