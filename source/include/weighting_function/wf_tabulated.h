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



#ifndef __WEIGHTING_FUNCTION_TABULATED_H
#define __WEIGHTING_FUNCTION_TABULATED_H

#include "weighting_function.h"


/* --- TabulatedWeightingFunctionWithWall --- */

template<typename T>
class MiniMatrix
{
 protected:
  size_t m_rows, m_columns;
  T *m_data;

  void alloc() {
    if (m_data)
      delete [] m_data;

    m_data = new T[m_rows*m_columns];
  }

 public:
  MiniMatrix(): m_data(NULL) {}
  ~MiniMatrix() {
    if (m_data)
      delete [] m_data;
  }

  void setSize(size_t rows, size_t columns) {
    m_rows = rows; m_columns = columns;
    alloc();
  }

  T &operator()(size_t y, size_t x) {
    return m_data[y*m_columns+x];
  }
  const T &operator()(size_t y, size_t x) const {
    return m_data[y*m_columns+x];
  }
};


class TabulatedWeightingFunctionWithWall: public WeightingFunctionWithWall
{
 protected:
  double m_rc_inv;
  WeightingFunction *m_wf;

  size_t m_n_bins;
  string m_wf_name;

  MiniMatrix<double> m_interpolate_table;
  MiniMatrix<double> m_weight_table;
  MiniMatrix<point_t> m_local_gradient_table;

  void init();

  virtual void setup();

 public:
  TabulatedWeightingFunctionWithWall(Node *parent);
  virtual ~TabulatedWeightingFunctionWithWall();

  virtual double interpolateWithDist(const Pairdist &r, const point_t &normal, double dist) const {
    return
      m_interpolate_table((size_t) (r.abs()*m_n_bins*m_rc_inv), (size_t) (dist*(m_n_bins-1)));
  }
  virtual double weightWithDist(const Pairdist &r, const point_t &normal, double dist) const {
    return
      m_weight_table((size_t) (r.abs()*m_n_bins*m_rc_inv), (size_t) (dist*(m_n_bins-1)));
  }
  virtual point_t localGradientWithDist(const Pairdist &r, const point_t &normal, double dist) const {
    return
      m_local_gradient_table((size_t) (r.abs()*m_n_bins*m_rc_inv), (size_t) (dist*(m_n_bins-1)));
  }
};


#endif
