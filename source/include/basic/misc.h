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



#ifndef __MISC_H
#define __MISC_H

#include "geometric_primitives.h"

typedef math_vector_t<int> int_point_t;
typedef math_vector_t<bool> bool_point_t;

#define MAX_COLOURS 10


/* Loop over a vector, list, etc. */

#define FOR_EACH(type, var, code)               \
{                                               \
  type::iterator __end = var.end();             \
  for (type::iterator __iFE = var.begin();          \
       __iFE != __end; ++__iFE) {                       \
    code                                        \
  }                                             \
}


#define FOR_EACH_CONST(type, var, code)         \
{                                               \
  type::const_iterator __end = var.end();       \
  for (type::const_iterator i = var.begin();    \
       i != __end; ++i) {                       \
    code                                        \
  }                                             \
}


#endif
