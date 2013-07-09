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



#ifndef __LINEAR_REGRESSION_H
#define __LINEAR_REGRESSION_H

#include "calc.h"
#include "postprocessor.h"

class LinearRegression: public Postprocessor
{
protected:
  DataFormat *m_input_format;
  list<data_sp> m_data;
  int m_idx_x, m_idx_y;

  fn m_prefn, m_postfn;
    
  string m_xcol, m_ycol, m_str_prefn, m_str_postfn;

  void init();

public:
  LinearRegression(Node *parent, Simulation *simulation);
  ~LinearRegression();

  virtual void setup();

  virtual void describeInput(DataFormat *input_format);

  virtual void push(data_sp data);

  virtual void flush();
};

#endif
