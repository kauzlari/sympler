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



#ifndef __FP_VECTOR_H
#define __FP_VECTOR_H

#include <string>

using namespace std;

#include "typed_value.h"

#include "geometric_primitives.h"

/*!
 * Base class for all vector quantities
 */
class FPVector: public TypedValue
{
 public:
  /*!
   * Constructor
   */
  FPVector(string name);

  /*!
   * Destructor
   */
  virtual ~FPVector();

  /*!
   * Return the value
   */
  virtual Variant value() const = 0;

  /*!
   * Return a C expression for this scalar
   */
  virtual Variant toC() const = 0;
};

// next is currently not used; if necessary, toC() should be corrected similarly to FPScalarConstant
#if 0
/*!
 * A constant vector
 */
class FPVectorConstant: public FPVector
{
 public:
  /*!
   * The value of this constant
   */
  point_t m_value;

 public:
  /*!
   * Constructor
   */
  FPVectorConstant(string name, point_t value);

  /*!
   * Destructor
   */
  virtual ~FPVectorConstant();

  /*!
   * Return the value
   */
  virtual Variant value() const {
    Variant v(Variant::VECTOR);

    for (int i = 0; i < SPACE_DIMS; ++i)
      v.vector(i) = m_value[i];


    return v;
  }

  /*!
   * Return a C expression for this constant
   */
  virtual Variant toC() const {
    Variant v(Variant::VECTOR_STRING);

    for (int i = 0; i < SPACE_DIMS; ++i)
      v.vectorString(i) = ObjToString(m_value[i]);


    return v;
  }
};
#endif

/*! 
 * A vector variable
 */
class FPVectorVariable: public FPVector
{
 protected:
  /*!
   * A pointer to the point_t holding the vector
   */
  point_t *m_value;

  /*!
   * The expression for C export
   */
  string m_c_expr[SPACE_DIMS];

 public:
  /*!
   * Constructor
   */
  FPVectorVariable(string name, point_t *value, string *c_expr = NULL);

  /*!
   * Destructor
   */
  virtual ~FPVectorVariable();

  /*!
   * Return the value
   */
  virtual Variant value() const {
    Variant v(Variant::VECTOR);

    for (int i = 0; i < SPACE_DIMS; ++i)
      v.vector(i) = (*m_value)[i];

    return v;
  }

  /*!
   * Return a C expression for this value
   */
  virtual Variant toC() const {
    Variant v(Variant::VECTOR_STRING);

    for (int i = 0; i < SPACE_DIMS; ++i)
      v.vectorString(i) = m_c_expr[i];

    return v;
  }
};


#endif
