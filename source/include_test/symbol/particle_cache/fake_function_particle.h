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



#ifndef __FAKE_FUNCTION_PARTICLE_H
#define __FAKE_FUNCTION_PARTICLE_H 

#include "function_particle.h"


/*!
 * Fake class for function particle. Current(2018-04-20) purpose:
 * - replace operator() by simple non-runtime-compiled function
 */
class FakeFunctionParticle: public FunctionParticle
{
 /* protected: */
                                        
 public:
  /*!
   * Constructor
   */
  FakeFunctionParticle();
  
  /*!
   * Destructor
   */
  virtual ~FakeFunctionParticle();

  /*!
   * Call the function
   * @param[out] val Return value
   * @param[in] p Particle the function takes as the argument
   */
  virtual void operator()(void* val, Particle* p) const {
    switch(m_returnType) {

    case Variant::SCALAR : 
      *((double*) val) =  3.;
      break;

    case Variant::VECTOR : 
      for(size_t i = 0; i < SPACE_DIMS; ++i) 
	(*((point_t*) val))[i] = i;
      break;

    case Variant::TENSOR :
      for(size_t i = 0; i < SPACE_DIMS; ++i) {
	for(size_t j = 0; j < SPACE_DIMS; ++j) {
	  (*((tensor_t*) val)) (i,j) = j + i*SPACE_DIMS;
	}
      }
      break;

    default:
      throw gError
          ("FakeFunctionParticle::operator()",
           "Don't know how to handle return type "
	   + ObjToString(m_returnType));

    }
  }
};

#endif
