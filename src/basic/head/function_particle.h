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



#ifndef __FUNCTION_PARTICLE_H
#define __FUNCTION_PARTICLE_H 

#include "function_arbitrary.h"
// #include "colour_pair.h"
// #include "weighting_function.h"
#include "particle.h"


/*!
 * A function which takes a particle as its argument
 */
class FunctionParticle: public FunctionArbitrary
{
 protected:
  /*!
   * The particle color
   */
  size_t m_colour;

//  /*!
//   * The \a ColourPair , the colour of the particles belongs to
  //   */
//  ColourPair* m_cp;                                    
                                        
public:
  /*!
   * Constructor
   */
  FunctionParticle();

  /*!
   * Destructor
   */
  virtual ~FunctionParticle();

  /*!
   * Compile to C code
   */
  virtual void compile();

  /*!
   * Set the color pair and the color
   */
  void setColour/*Pair*/(/*ColourPair* cp, */size_t colour) {
//     m_cp = cp;
    m_colour = colour;
  }

  /*!
   * Call the function
   * @param val Return value
   * @param p Particle the function takes as the argument
   */
  void operator()(void* val, Particle* p) const {
    m_compiler.fn()(val, p, p->tag.data());
  }
};

#endif
