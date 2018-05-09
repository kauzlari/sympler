/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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



#ifndef __PARTICLE_VELS_H
#define __PARTICLE_VELS_H

#include "particle_cache_arbitrary.h"

#include "function_particle.h"
#include "symbol.h"


/*!
 * Will not introduce any new Symbol but set the velocity of each 
 * particle to the value computed by the given expression.
 */
class ParticleVels : public ParticleCacheArbitrary
{
 protected:
  
  /*!
   * Initialise the PropertyList.
   */
  void init();
  
  /*!
   * Helper for setting m_offset. This one overwrites 
   * \a ParticleCacheArbitrary::setupOffset, since 
   * \a ParticleVels modifies \a Particle::v and hence does not 
   * require any offset, hence no usage of 
   * \a Particle::s_tag_format[m_colour] methods such as 
   * addAttribute(..)
   */
  virtual void setupOffset();
  
  /*!
   * Helper function for polymorphic copying
   */
  virtual ParticleCache* copyMySelf()
  {
    return new ParticleVels(*this);
  }

  /*!
   * Helper function for setting the return type of \a m_function
   */
  virtual void setFunctionReturnType(){
    m_function->setReturnType(Variant::VECTOR);
  }
  
 public:

  /*!
   * Constructor
   */
  ParticleVels(/*Node*/Simulation* parent);
  
  /*!
   * Destructor
   */
  virtual ~ParticleVels();
  
  /*!
   * Setup this ParticleCache
   */
  virtual void setup();
  
  /*!
   * Compute the cache for particle \a p
   * @param[out] p The particle to compute a new velocity for
   */
  virtual void computeCacheFor(Particle* p) {
    
    (*m_function)(&(p->v), p);            
  }
  
  /*!
   * Register the additional degrees of freedom with the 
   * \a Particle s \a DataFormat
   * FIXME: Somehow, most of the classes do not use this function 
   * anymore. So CHECK if it is:
   * - A: Deprecated?
   * - B: A useful refactoring that you forgot to use? (Note that 
   * there is also the more recently created 
   * ParticleCacheArbitrary::setupOffset()!)
   */
  virtual void registerWithParticle()
  {
  }
  
  /*!
   * Is this \a ParticleCache identical to the given one?
   * @param c \a ParticleCache to compare to
   * FIXME: Essentially copy&pasted from \a ParticleVector, so is
   * there some refactoring possible? Why does the parent 
   * \a ParticleCacheArbitrary not define anything? Because of some 
   * nasty children?
   */
  virtual bool operator==(const ParticleCache &c) const
  {
    if (typeid(c) == typeid(*this)) 
      {
	ParticleVels *cc = (ParticleVels*) &c;
	
	return
	(
	 m_expression == cc->m_expression &&
	 m_overwrite == cc->m_overwrite &&
	 m_stage == cc->m_stage &&
	 m_colour == cc->m_colour &&
	 m_symbolName == cc->m_symbolName &&
	 m_datatype == cc->m_datatype &&
	 m_offset == cc->m_offset
	 );
      }
    else return false;
  }
  
};

#endif
