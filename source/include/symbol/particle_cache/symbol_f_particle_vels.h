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


#ifndef __SYMBOL_F_PARTICLE_VELS_H
#define __SYMBOL_F_PARTICLE_VELS_H

#include "symbol_f_particle_arbitrary.h"

using namespace std;

class Simulation;

/*!
 * Per particle translational force called as a \a Symbol
 */
class SymbolFParticleVels : public SymbolFParticleArbitrary
{

 protected:  

  /*!
   * Initialise the PropertyList.
   */
  void init();

  /*!
   * Helper for setting m_offset. This one uses 
   * \a SymbolFParticleArbitrary::setupOffset, which in turn overrides 
   * \a ParticleCacheArbitrary::setupOffset, since 
   * \a SymbolFParticleVels modifies \a Particle::force and hence does 
   * not require any offset, hence no usage of 
   * \a Particle::s_tag_format[m_colour] methods such as 
   * addAttribute(..)
   */
  virtual void setupOffset();

  /*!
   * Helper function for polymorphic copying
   */
  virtual ParticleCache* copyMySelf()
  {
    return new SymbolFParticleVels(*this);
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
   * @param simulation Pointer to the \a Simulation object
   */
  SymbolFParticleVels(Simulation *simulation);
  
  /*!
   * Destructor
   */
  virtual ~SymbolFParticleVels();
  
  /*!
   * Setup this \a Symbol
   */
  virtual void setup();
  
  /*!
   * Compute the cache for particle \a p
   * @param[out] p The particle to compute a new force for
   */
  virtual void computeCacheFor(Particle* p) {
    
    (*m_function)(&(p->force[m_offset]), p);            
  }
        
};

#endif
