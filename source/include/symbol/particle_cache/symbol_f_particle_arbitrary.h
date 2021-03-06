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


#ifndef __SYMBOL_F_PARTICLE_ARBITRARY_H
#define __SYMBOL_F_PARTICLE_ARBITRARY_H

#include "particle_cache_arbitrary.h"

#include "function_particle.h"
#include "simulation.h"

using namespace std;

class Simulation;

/*!
 * Per particle force/flux called as a \a Symbol
 */
class SymbolFParticleArbitrary : public ParticleCacheArbitrary
{

 protected:  

  /*!
   * Initialise the PropertyList.
   */
  void init();
  
  /*!
   * Helper for setting m_offset. This one overwrites 
   * \a ParticleCacheArbitrary::setupOffset, since any 
   * \a SymbolFParticleArbitrary modifies an existing force/flux and 
   * hence does, e.g., not add symbols with 
   * \a Particle::s_tag_format[m_colour].addAttribute(..)
   */
  virtual void setupOffset();
  
  /*!
   * Helper function for polymorphic copying
   */
  virtual ParticleCache* copyMySelf() = 0;

  virtual void checkIfSymbolNameUntouched() const;
 
 public:
  
  /*!
   * Constructor
   * @param simulation Pointer to the \a Simulation object
   */
  SymbolFParticleArbitrary(Simulation *simulation);
  
  /*!
   * Destructor
   */
  virtual ~SymbolFParticleArbitrary();
  
  /*!
   * Setup this \a Symbol
   */
  virtual void setup();

  /*!
   * Possibility for the Node to precompute stuff before \a Symbol s 
   * and \a Force s start computing
   */
  virtual void precompute()
  {
    m_offset = ((Simulation*) m_parent) -> controller() -> forceIndex(); 
  }
  
  /*!
   * Compute the cache for particle \a p
   * @param[out] p The particle to compute a new force/flux for
   */
  virtual void computeCacheFor(Particle* p) = 0;
  
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
  virtual void registerWithParticle() {}
  
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
	SymbolFParticleArbitrary *cc = (SymbolFParticleArbitrary*) &c;
	
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

  /*!
   * Currently (2018-05-25) this function is only used for unittesting
   */
  virtual void setForceIndexTo(size_t forceIndex) {
    m_offset = forceIndex;
  }
  
};

#endif
