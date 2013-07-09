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



#ifndef __TYPED_VALUE_H
#define __TYPED_VALUE_H

#include <string>

using namespace std;

#include "function_node.h"

/*!
 * Represents a constant or a variable within the
 * parse tree
 */
class TypedValue: public FunctionNode
{
 public:
  enum typed_value_t {
    SCALAR = 0,
    VECTOR = 1,
    TENSOR = 2,
    EO_TYPED_VALUE = 3
  };

 protected:
  /*!
   * The type of this specific node
   */
  typed_value_t m_type;

 public:
  /*!
   * Constructor
   */
  TypedValue(string name, typed_value_t t);

  /*!
   * Destructor
   */
  virtual ~TypedValue();
};


#endif
