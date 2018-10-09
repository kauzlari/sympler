/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018,
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#ifndef __MATH_HELPER_H
#define __MATH_HELPER_H

class MathHelper {

	public:

		/*!
		 * Computes matrix-matrix product result_ab = mat1_ac . mat2_cb where
		 * result is size_a x size_b, mat1 is size_a x size_c, and mat2 is
		 * size_c x size_b. Entries of result do not have to be initialised to 0.
		 * Memory for mat1, mat2, result must have been allocated. Ordering is
		 * row-major.
		 * @param[in] size_a Number of rows of mat1 and result
		 * @param[in] size_b Number of columns of mat2 and result
		 * @param[in] size_c Number of columns of mat1 and rows of mat2
		 * @param[in] mat1 Handle to row-major storage of mat1
		 * @param[in] mat2 Handle to row-major storage of mat2
		 * @param[out] result Handle to row-major storage of result
		 */
		static void matrixMatrixProd
		(size_t size_a, size_t size_b, size_t size_c, const double* mat1,
				const double* mat2, double* result) {

			for (size_t a = 0; a < size_a; ++a) {

				for (size_t b = 0; b < size_b; ++b) {

					size_t idx = b + size_b * a;
					result[idx] = 0.;

					for (size_t c = 0; c < size_c; ++c) {

						result[idx] += mat1[c + size_c * a] * mat2[b + size_b * c];

					} // end of for (size_t c...

				} // end of for (size_t b...

			} // end of for (size_t a...

		}

		/*!
		 * Computes matrix-matrix product result_ab = mat1_ac . mat2_cb^T, i.e.,
		 * i.e., mat2 is transposed before multiplication.
		 * result is size_a x size_b, mat1 is size_a x size_c, and mat2 is
		 * size_c x size_b. Entries of result do not have to be initialised to 0.
		 * Memory for mat1, mat2, result must have been allocated. Ordering is
		 * row-major.
		 * @param[in] size_a Number of rows of mat1 and result
		 * @param[in] size_b Number of columns of mat2 and result
		 * @param[in] size_c Number of columns of mat1 and rows of mat2
		 * @param[in] mat1 Handle to row-major storage of mat1
		 * @param[in] mat2 Handle to row-major storage of mat2
		 * @param[out] result Handle to row-major storage of result
		 */
		static void matrixMatrixTProd
		(size_t size_a, size_t size_b, size_t size_c, const double* mat1,
				const double* mat2, double* result) {

			for (size_t a = 0; a < size_a; ++a) {

				for (size_t b = 0; b < size_b; ++b) {

					size_t idx = b + size_b * a;
					result[idx] = 0.;

					for (size_t c = 0; c < size_c; ++c) {

						result[idx]
									 += mat1[c + size_c * a]
													 * mat2[c + size_c * b]; // transposition!

					} // end of for (size_t c...

				} // end of for (size_t b...

			} // end of for (size_t a...
		}

};

#endif
