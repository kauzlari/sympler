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



#ifndef __OUTPUT_DCD_H
#define __OUTPUT_DCD_H

#include <fstream>

#include "output.h"

using namespace std;


/*!
 * Output module for DCD files. This has only been tested with VMD 
 * (http://www.ks.uiuc.edu/Research/vmd/) right now and cannot
 * be used to export degrees of freedom.
 */
class OutputDCD: public Output
{
protected:
  /*!
   * Format description of the data passed on by the parent
   * \a Postprocessor
   */
  DataFormat *m_input_format;
	
  /*!
   * File name
   */
  string m_fn;
	
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
   * The number of atoms in the simulation. Note: This is initialized from
   * the first time step. Changing particle numbers are not supported by
   * the DCD format.
   */
  size_t m_n_atoms;

  /*!
   * The output stream
   */
  ofstream m_s;
	
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
  OutputDCD(Node *parent, Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~OutputDCD();

  /*!
   * Get the index of the "__data_format" field
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

  /*!
   * Closes the stream
   */
  virtual void flush();
};

#endif
