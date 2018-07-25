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



#ifndef __FP_TENSOR_H
#define __FP_TENSOR_H

#include <string>

using namespace std;

#include "typed_value.h"

#include "geometric_primitives.h"

/*!
 * Base class for all tensor quantities
 */
class FPTensor: public TypedValue
{
  public:
  /*!
   * Constructor
   */
    FPTensor(string name);

  /*!
     * Destructor
   */
    virtual ~FPTensor();

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
 * A constant tensor
 */
class FPTensorConstant: public FPTensor
{
  public:
  /*!
   * The value of this constant
   */
    tensor_t m_value;

  public:
  /*!
   * Constructor
   */
    FPTensorConstant(string name, tensor_t value);

  /*!
     * Destructor
   */
    virtual ~FPTensorConstant();

  /*!
     * Return the value
   */
    virtual Variant value() const {
      Variant v(Variant::TENSOR);

      for (int i = 0; i < SPACE_DIMS; ++i)
        for (int j = 0; j < SPACE_DIMS; ++j)
          v.tensor(i, j) = m_value(i, j);

      return v;
    }

  /*!
     * Return a C expression for this constant
   */
    virtual Variant toC() const {
      Variant v(Variant::TENSOR_STRING);
      for (int i = 0; i < SPACE_DIMS; ++i)
        for (int j = 0; j < SPACE_DIMS; ++j)
          v.tensorString(i, j) = ObjToString(m_value(i, j));

      return v;
    }
};

#endif

/*! 
 * A tensor variable
 */
class FPTensorVariable: public FPTensor
{
  protected:
  /*!
   * A pointer to the tensor_t holding the tensor
   */
    tensor_t *m_value;

  /*!
     * The expression for C export
   */
    string m_c_expr[SPACE_DIMS][SPACE_DIMS];

  public:
  /*!
   * Constructor
   */
    FPTensorVariable(string name, tensor_t *value, string c_expr[SPACE_DIMS][SPACE_DIMS] = NULL);

  /*!
     * Destructor
   */
    virtual ~FPTensorVariable();

  /*!
     * Return the value
   */
    virtual Variant value() const {
      Variant v(Variant::TENSOR);

      for (int i = 0; i < SPACE_DIMS; ++i)
        for (int j = 0; j < SPACE_DIMS; ++j)
          v.tensor(i, j) = (*m_value)(i, j);

      return v;
    }

  /*!
     * Return a C expression for this value
   */
    virtual Variant toC() const {
      Variant v(Variant::TENSOR_STRING);

      for (int i = 0; i < SPACE_DIMS; ++i)
        for (int j = 0; j < SPACE_DIMS; ++j)
          v.tensorString(i, j) = m_c_expr[i][j];

      return v;
    }
};

#endif
