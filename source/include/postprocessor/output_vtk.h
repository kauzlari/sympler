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



#ifndef __OUTPUT_VTK_H
#define __OUTPUT_VTK_H

#include <fstream>

#include "output.h"

using namespace std;


/*!
 * Write the data to a VTK file (Visualization ToolKit). This module
 * can handle particle positions and \a GridAverager output and is
 * the preferred method for writing data to disk as VTKs
 * are suffiently
 */
class OutputVTK: public Output
{
protected:
  /*!
   * Format description of the data passed on by the parent
   * \a Postprocessor
   */
  DataFormat *m_input_format;
	
  /*!
   * The name of the VTK to write
   */
  string m_fn;
  
  /*!
   * Format (either 'binary' or 'ascii')
   */
  string m_format;
	
  /*!
   * The output step
   */
  int m_step_counter;

  /*!
   * The index of the "__data_format" entry in the data passed on from the parent
   * \a Postprocessor
   */
  int m_idx_data_format;

  /*!
   * The index that contains the number of cells (only valid if data is passed
   * from a \a GridAverager)
   */
  int m_idx_n_cells;

  /*!
   * Initialize the property list
   */
  void init();
	
public:
  /*!
   * Constructor
   * @param parent Parent node in the \a Postprocessor chain
   * @param simulation Pointer to the main \a Simulation object
   */
  OutputVTK(Node *parent, Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~OutputVTK();

  /*!
   * Get the index of the "__data_format" field and check whether we are
   * getting particle positions or a grid
   */ 
  virtual void setup();
	
  /*!
   * Called by the parent to pass on its \a DataFormat
   */
  virtual void describeInput(DataFormat *input_format);
	
  /*!
   * Data from the parent \a Postprocessor or \a Meter
   */
  virtual void push(data_sp data);
};

#endif
