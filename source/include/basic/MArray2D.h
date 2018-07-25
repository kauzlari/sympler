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



#ifndef __DENSEMAT_MAR_H_
#define __DENSEMAT_MAR_H_

#include <general.h>
#ifdef WITH_ARRAY_TYPES

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <ctype.h>
#include "tnt/tnt.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <vector>

#define MMBanner "%%MatrixMarket"
#define MM_MAX_LINE_LENGTH 1025
#define MM_MAX_TOKEN_LENGTH 64
#define MM_DENSE_STR  "array"
#define MM_REAL_STR	  "real"
#define MM_MTX_STR	  "matrix"
#define MM_GENERAL_STR  "general"

TNT::Array2D< double>* read_MMheader(string path);

/*!Class for importing and handling matrices in MatrixMarket format*/
class MArray2D {

	std::string path;
public:
	TNT::Array2D<double> *p_array2d;
	MArray2D();
	MArray2D(std::string filename);
	MArray2D(const MArray2D &o);
	virtual ~MArray2D();
	MArray2D(int r, int c);
	MArray2D(int r, int c, int z);
	size_t dim1() {
		return p_array2d ? p_array2d->dim1() : 0;
	}
	size_t dim2() {
		return p_array2d ? p_array2d->dim2() : 0;
	}
	MArray2D scalarmult(double s);
	MArray2D& row(int i);
	MArray2D& column(int j);
	inline double operator()(size_t i, size_t j) const{
		return (*p_array2d)[i][j];
	}
	inline double& operator()(size_t i, size_t j) {
			return (*p_array2d)[i][j];
		}
	std::string toMTXString(const std::string &eol = "\n");
#ifdef HAVE_JAMA_JAMA_LU_H
	friend MArray2D &invert(const MArray2D &m1);
#endif
	friend MArray2D row(MArray2D &m, int i);
	friend MArray2D column(MArray2D &m, int j);
	inline double* operator[](int i) {
		return (*p_array2d)[i];
	}



	/*!
	 * Assigment operator
	 * @param copy_marray2d Copy from this object
	 */
	MArray2D &operator=(const MArray2D &copy_marray2d);
};
std::ostream& operator<<(std::ostream &s, MArray2D &m);
MArray2D operator+(const MArray2D &m1, const MArray2D &m2);
MArray2D operator-(const MArray2D &m1, const MArray2D &m2);
MArray2D matmult(const MArray2D &m1, const MArray2D &m2);

#endif /*WITH_ARRAY_TYPES*/
#endif /*DENSEMAT_MAR_H_*/
