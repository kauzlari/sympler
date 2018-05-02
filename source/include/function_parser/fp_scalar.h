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



#ifndef __FP_SCALAR_H
#define __FP_SCALAR_H

#include <string>

using namespace std;

#include "typed_value.h"

/*!
 * Base class for all scalar quantities
 */
class FPScalar: public TypedValue
{
 public:
  /*!
   * Constructor
   */
  FPScalar(string name);

  /*!
   * Destructor
   */
  virtual ~FPScalar();
};


/*!
 * A scalar constant
 */
class FPScalarConstant: public FPScalar
{
 public:
  /*!
   * The value of this constant
   */
  double m_value;

 public:
  /*!
   * Constructor
   */
  FPScalarConstant(string name, double value);

  /*!
   * Destructor
   */
  virtual ~FPScalarConstant();

  /*!
   * Return the value
   */
  virtual Variant value() const {
    Variant v(Variant::SCALAR);

    v.scalar() = m_value;

    return v;
  }

  /*!
   * Return a C expression for this constant
   */
  virtual Variant toC() const {
    Variant v(Variant::SCALAR_STRING);

    // lets be sure that each integer value is compiled as a double
    if(m_value - int(m_value) == 0)
      v.scalarString() = "(" + ObjToString(int(m_value)) + ".0)";
    else
      v.scalarString() = m_name /* OLD STYLE can produce undesired integer divisions like "1/8": ObjToString(m_value)*/;

    return v;
  }
};


/*! 
 * A scalar variable
 */
class FPScalarVariable: public FPScalar
{
 protected:
  /*!
   * A pointer to the double holding the value
   */
  double *m_value;

  /*!
   * The expression for C export
   */
  string m_c_expr;

 public:
  /*!
   * Constructor
   */
  FPScalarVariable(string name, double *value, string c_expr = "");

  /*!
   * Destructor
   */
  virtual ~FPScalarVariable();

  /*!
   * Return the value
   */
  virtual Variant value() const {

    Variant v(Variant::SCALAR);

    // FIXME: currently (2018-02-12) this check is only done here, but
    // should be generally useful to avoid seg-faults. Put also at same
    // places for vectors and tensors, and briefly think if there is a
    // generalising refactoring that could be done to avoid code-
    // repetition.
    if(!m_value)
      throw gError("FPScalarVariable::value()", "ERROR: Tried to return value of my non-initialised double* m_value = NULL.");
    
    v.scalar() = *m_value;

    return v;
  }

  /*!
   * Return a C expression for this constant
   */
  virtual Variant toC() const {
    Variant v(Variant::SCALAR_STRING);

    v.scalarString() = m_c_expr;

    return v;
  }
};


#endif
