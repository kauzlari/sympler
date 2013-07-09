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



#include "fp_tensor.h"


/* FPTensor */

FPTensor::FPTensor(string name): TypedValue(name, TENSOR)
{
}


FPTensor::~FPTensor()
{
}

// commented out because currently not used
#if 0
/* FPTensorConstant */

FPTensorConstant::FPTensorConstant(string name, tensor_t value)
  : FPTensor(name), m_value(value)
{
}


FPTensorConstant::~FPTensorConstant()
{
}
#endif

/* FPTensorVariable */

FPTensorVariable::FPTensorVariable(string name, tensor_t *value, string c_expr[SPACE_DIMS][SPACE_DIMS])
  : FPTensor(name), m_value(value)
{
  if (!c_expr) {
//     MSG_DEBUG("FPTensorVariable::FPTensorVariable", "!c_expr-case");
    for (int i = 0; i < SPACE_DIMS; ++i)
      for (int j = 0; j < SPACE_DIMS; ++j)
        m_c_expr[i][j] = name + ObjToString(i) + "_" + ObjToString(j);
  } 
  else {
    for (int i = 0; i < SPACE_DIMS; ++i)
    {
      for (int j = 0; j < SPACE_DIMS; ++j)
      {
        m_c_expr[i][j] = c_expr[i][j];
//         MSG_DEBUG("FPTensorVariable::FPTensorVariable", "c_expr-case: c_expr[" << i << "][" << j << "] = " << c_expr[i][j]);
      }
    }
  }
}


FPTensorVariable::~FPTensorVariable()
{
}

