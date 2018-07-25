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



#include "fp_vector.h"


/* FPVector */

FPVector::FPVector(string name): TypedValue(name, VECTOR)
{
}


FPVector::~FPVector()
{
}

// commented out because currently not used
#if 0
/* FPVectorConstant */

FPVectorConstant::FPVectorConstant(string name, point_t value)
  : FPVector(name), m_value(value)
{
}


FPVectorConstant::~FPVectorConstant()
{
}
#endif

/* FPVectorVariable */

FPVectorVariable::FPVectorVariable(string name, point_t *value, string *c_expr)
  : FPVector(name), m_value(value)
{
  if (!c_expr) {
    for (int i = 0; i < SPACE_DIMS; ++i)
      m_c_expr[i] = name + ObjToString(i);
  } else {
    for (int i = 0; i < SPACE_DIMS; ++i)
      m_c_expr[i] = c_expr[i];
  }
}


FPVectorVariable::~FPVectorVariable()
{
}

