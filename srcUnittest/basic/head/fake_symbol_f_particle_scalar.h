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


#ifndef __FAKE_SYMBOL_F_PARTICLE_SCALAR_H
#define __FAKE_SYMBOL_F_PARTICLE_SCALAR_H

#include "symbol_f_particle_scalar.h"

using namespace std;

/*!
 * Fake class for \a SymbolFParticleScalar. Current purpose(s):
 * - public version of protected method 
 *   SymbolFParticleScalar::setupOffset()
 */
class FakeSymbolFParticleScalar : public SymbolFParticleScalar
{

 public:
  
  /*!
   * Constructor
   * @param simulation Pointer to the \a Simulation object
   */
  FakeSymbolFParticleScalar(Simulation *simulation);
  
  /*!
   * Destructor
   */
  virtual ~FakeSymbolFParticleScalar();

  /*!
   * Performs the same task as the protected method 
   * \a SymbolFParticleScalar::setupOffset() but without setting 
   * \a SymbolFParticleArbitrary::m_offset as well, hence without 
   * calling \a SymbolFParticleArbitrary::setupOffset(). The attribute 
   * added to \a Particle::s_tag_format is special and hard-coded (see 
   * function definition) for test purposes only.
   */
  void setupForceOffsetForTestsOnly();

};

#endif
