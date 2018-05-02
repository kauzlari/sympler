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



#ifndef __BINARY_OPERATORS
#define __BINARY_OPERATORS

#include "function_node.h"


/*!
 * FNBinaryOperator is the base class for all binary operators
 */
class FNBinaryOperator: public FunctionNode
{
 protected:					
  /*!
   * The two parameters of the binary operator
   */
  FunctionNode *m_a, *m_b;

 public:
  /*!
   * Constructor
   */
  FNBinaryOperator(string name);

  /*!
   * Destructor
   */
  virtual ~FNBinaryOperator();

  /*!
   * Set the value of the two parameters
   */
  virtual void setBinary(FunctionNode *a, FunctionNode *b);
};



/* MACROS */

#define FOR_EACH_DOUBLE2(where, varianta, variantb, code)		\
  {									\
    Variant::variant_type_t t;						\
    Variant vca = varianta;						\
    Variant vcb = variantb;						\
									\
    t = vca.sameType(variantb);						\
									\
    if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) {\
      Variant r(t);							\
      vector_double_t &vr = r.doubles();				\
      vector_double_t va = vca.doubles();				\
      vector_double_t vb = vcb.doubles();				\
									\
      for (size_t i = 0; i < va.size(); ++i) {				\
	code;								\
      }									\
									\
      return r;								\
    } else								\
      throw gError							\
	(#where,							\
	 "Can't handle Variant type id " + ObjToString(t) + ".");	\
  }


#define FOR_EACH_STRING2(where, varianta, variantb, code)		\
  {									\
    Variant::variant_type_t t;						\
    Variant vca = varianta;						\
    Variant vcb = variantb;						\
									\
    t = vca.sameType(variantb);						\
									\
    if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	\
	t == Variant::TENSOR_STRING) {					\
      Variant r(t);							\
      vector_string_t &vr = r.strings();				\
      vector_string_t va = vca.strings();				\
      vector_string_t vb = vcb.strings();				\
									\
      for (size_t i = 0; i < va.size(); ++i) {				\
	code;								\
      }									\
									\
      return r;								\
    } else								\
      throw gError							\
	(#where,							\
	 "Can't handle Variant type id " + ObjToString(t) + ".");	\
  }


#endif
