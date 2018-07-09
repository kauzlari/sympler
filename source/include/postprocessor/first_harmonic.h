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



#ifndef __FIRST_HARMONIC_H
#define __FIRST_HARMONIC_H

#include "postprocessor.h"

/*!
 * Calculates the first harmonic (i.e., the amplitude of the first Fourier mode).
 * The input data for this \a Postprocessor needs three entries:
 *   - The current time
 *   - A column containing the x-positions of the data points
 *   - A column containing the f(x)-information to calculate the first harmonic for
 */
class FirstHarmonic: public Postprocessor
{
protected:
  /*!
   * Format description of the data passed on by the parent
   * \a Postprocessor
   */
  DataFormat *m_input_format;

  /*!
   * The index of the time information
   */
  int m_idx_time;

  /*!
   * The index of the column containing the x-positions
   */
  int m_idx_x;

  /*!
   * The index of the column containing the f(x)-information
   */
  int m_idx_y;

  /*!
  int m_idx_left;
  int m_idx_right;
  int m_idx_grid_size;
  */

  // next two used for abort, when the noise becomes stronger than 
  // the exponential decay

  /*!
   *
   */
  double oldValue;

  /*!
   *
   */
  bool negativeFound;

  void init();
  string m_species;
      
public:
  FirstHarmonic(Node *parent, Simulation *simulation);
  ~FirstHarmonic();

  virtual void setup();

  virtual void describeInput(DataFormat *input_format);

  virtual void push(data_sp data);    
};

#endif
