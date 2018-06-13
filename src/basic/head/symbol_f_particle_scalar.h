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


#ifndef __SYMBOL_F_PARTICLE_SCALAR_H
#define __SYMBOL_F_PARTICLE_SCALAR_H

#include "symbol_f_particle_arbitrary.h"

using namespace std;

class Simulation;

/*!
 * Implementation of a \a Symbol modifying a force/flux on a scalar 
 * degree of freedom (introduced by an \a Integrator).
 * NOTE the re-interpretation of 
 * \a ParticleCacheArbitrary::m_symbolName: The name of the scalar 
 * from which we modify its force/flux
 * \a ParticleCacheArbitrary::m_offset: The offset of the force/flux on 
 * the scalar
 */
class SymbolFParticleScalar : public SymbolFParticleArbitrary
{

 protected:  

  /*!
   * Internal helper which is set once to the respective force indices 
   * of the force on \a m_variableName.
   * FIXME: when also implementing classes SymbolFParticleVector/Tensor,
   * this common member (and others?) justifies a common parent class.
   */
  size_t m_forceOffset[FORCE_HIST_SIZE];

  /*!
   * Name of the variable from which the force is modified. An 
   * additional string variable is used in order to keep the original 
   * meaning of \a Symbol::m_symbolName sharp.
   * FIXME: when also implementing classes SymbolFParticleVector/Tensor,
   * this common member (and others?) justifies a common parent class.
   */
  string m_variableName;
  
  /*!
   * Initialise the PropertyList.
   */
  void init();

  /*!
   * Helper for setting m_offset. This one uses 
   * \a SymbolFParticleArbitrary::setupOffset, which in turn overrides 
   * \a ParticleCacheArbitrary::setupOffset, since 
   * \a SymbolFParticleScalar modifies the force on the scalar and not 
   * the scalar itself. Therefore, the internal array size_t 
   * \a m_force_offset[FORCE_HIST_SIZE] is set to the respective force 
   * indices. This is then used to set \a m_offset dynamically to the 
   * right one in the array.
   */
  virtual void setupOffset();

  /*!
   * Helper function for polymorphic copying
   */
  virtual ParticleCache* copyMySelf()
  {
    return new SymbolFParticleScalar(*this);
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
   * @param simulation Pointer to the \a Simulation object
   */
  SymbolFParticleScalar(Simulation *simulation);
  
  /*!
   * Destructor
   */
  virtual ~SymbolFParticleScalar();
  
  /*!
   * Setup this \a Symbol
   */
  virtual void setup();
  
  /*!
   * Compute the flux/force for particle \a p
   * Note that the force index is stored in \a m_offset
   * @param[out] p The particle to compute a new flux for
   */
  virtual void computeCacheFor(Particle* p) {

    (*m_function)
      (&(p -> tag.doubleByOffset(m_forceOffset[m_offset])), p);

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
	SymbolFParticleScalar *cc = (SymbolFParticleScalar*) &c;

	bool tmpBool = true;
	for (size_t i = 0; i < FORCE_HIST_SIZE; ++i) {
	  tmpBool = tmpBool && (m_forceOffset[i] == cc->m_forceOffset[i]);
	  if (!tmpBool) break;
	}

	return(tmpBool &&
	       SymbolFParticleArbitrary::operator==(c) &&
	       (m_variableName == cc->m_variableName)
	       );
      }
    else return false;
  }

  /*!
   * Return \a m_variableName
   * @param[return] The string stored in \a m_variableName
   */
  virtual const string& returnVariableName() const {
    return m_variableName;
  }

  /*!
   * Return the specified forceOffset from \a m_forceOffset
   * @param[in] index The index of the requested forceOffset 
   * @param[return] The requested forceOffset
   */
  virtual size_t returnForceOffset(size_t index) const {
    return m_forceOffset[index];
  }
  
};

#endif
