/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


#ifndef __GSL_HELPER_H
#define __GSL_HELPER_H

#include <gsl/gsl_linalg.h>

// the following definitions are done because I do not know how to turn range-checking 
// off for the GSL
#ifndef GSL_VECTOR_SET
#define GSL_VECTOR_SET(v, i, x) (v)->data[i*(v)->stride] = x
#endif
#ifndef GSL_VECTOR_GET
#define GSL_VECTOR_GET(v, i) (v)->data[i*(v)->stride]
#endif
#ifndef GSL_MATRIX_SET
#define GSL_MATRIX_SET(m, i, j, x) (m)->data[i * (m)->tda + j] = x
#endif
#ifndef GSL_MATRIX_GET
#define GSL_MATRIX_GET(m, i, j) (m)->data[i * (m)->tda + j]
#endif


#endif
