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


#ifndef __PARTICLE_SCALAR_H
#define __PARTICLE_SCALAR_H

#include "particle_cache_arbitrary.h"

// #define PCA_MAX_STAGE 2

/*!
 * User-defined scalar symbol for a particle that depends only
 * on previously computed particle properties.
 */
class ParticleScalar : public ParticleCacheArbitrary
{
  protected:
  
    /*!
     * Initialise the PropertyList.
     */
    void init();
  
//    /*!
//     * Helper for setting m_offset
    //     */
//    void setupOffset();
    
    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new ParticleScalar(*this);
    }
    
    /*!
     * Helper function for setting the return type of \a m_function
     */
    virtual void setFunctionReturnType(){
      m_function->setReturnType(Variant::SCALAR);
    }

 public:
  /*!
   * Constructor
   */
    ParticleScalar(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ParticleScalar();

      /*!
     * Compute the cache for particle \a p
     * @param p The particle to compute values for
       */
    virtual void computeCacheFor(Particle* p)
    {
      
      (*m_function)(&(p->tag.doubleByOffset(m_offset)), p);
      
    }

  /*!
     * Register the additional degrees of freedom with the \a Particle s \a DataFormat
   */
    virtual void registerWithParticle()
    {
    }
    
  /*!
     * Are those two caches identical?
     * @param c Is this the same cache?
   */
    virtual bool operator==(const ParticleCache &c) const
    {
      if (typeid(c) == typeid(*this)) 
      {
        ParticleScalar *cc = (ParticleScalar*) &c;

        return
          (
          m_expression == cc->m_expression &&
          m_overwrite == cc->m_overwrite &&
          m_stage == cc->m_stage &&
          m_colour == cc->m_colour &&
          m_symbolName == cc->m_symbolName &&
          m_datatype == cc->m_datatype
          );
      }
      else return false;
    }
    
  /*!
     * Setup this ParticleCache
   */
    virtual void setup();
};

#endif
