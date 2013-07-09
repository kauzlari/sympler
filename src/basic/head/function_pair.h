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



#ifndef __FUNCTION_PAIR_H
#define __FUNCTION_PAIR_H 

#include "function_arbitrary.h"
#include "pairdist.h"

class ColourPair;


/*!
 * A function which takes a pair as its argument
 */
class FunctionPair: public FunctionArbitrary
{
protected:
  /*!
   * The \a ColourPair corresponding to the pairs of this function
   */
  ColourPair* m_cp;

public:
  /*!
   * Constructor
   */
  FunctionPair();

  /*!
   * Destructor
   */
  virtual ~FunctionPair();

  /*!
   * Compile the function
   */
  virtual void compile();

  /*!
   * Set the \a ColourPair
   * @param cp The \a ColourPair to use
   */
  void setColourPair(ColourPair* cp) {
    assert(cp);
    m_cp = cp;
  }

  /*!
   * Call the function
   * @param val Return value
   * @param pr Pairdist this function takes as the argument
   */
  void operator()(void* val, Pairdist* pr) const {
    m_compiler.fn()
      (val,
       &(pr->m_distance),
       pr->tag.data(),
       pr->firstPart(),
       pr->secondPart(),
       pr->firstPart()->tag.data(),
       pr->secondPart()->tag.data());
  }
};

#endif
