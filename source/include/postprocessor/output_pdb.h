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



#ifndef __OUTPUT_POS_VEL_PDB_H
#define __OUTPUT_POS_VEL_PDB_H

#include <fstream>

#include "output.h"

using namespace std;


/*!
 * Write the data to a PDB (Protein DataBase) file. Note: This type of file cannot store
 * velocities.
 */
class OutputPDB: public Output
{
 protected:
  /*!
   * Format description of the data passed on by the parent
   * \a Postprocessor
   */
  DataFormat *m_input_format;
	
  /*!
   * Name of the file write to
   */
  string m_posfn;
	
  /*!
   * Stream to write to
   */
  ofstream m_pos_s;
	
  /*!
   * Use a single file for every time step. Fixme!!! Does it make sense to use
   * just a large file.
   */
  bool m_multiple_files;
	
  /*!
   * The output step
   */
  int m_step_counter;

  /*!
   * The index where to find the positions in the data passed
   * to this output module.
   */
  int m_idx_positions;
	
  /*!
   * Write a list of vectors to the output file
   */
  void writeVectorPoint(ostream &s, vector_point_sp vp);
	
  /*!
   * Initialize the property list
   */
  void init();
	
  /*!
   * Convert a \a point_t to a \a string
   */
  static string point2string(const point_t &p, bool comma = true);
	
 public:
  /*!
   * Constructor
   * @param parent Parent node in the \a Postprocessor chain
   * @param simulation Pointer to the main \a Simulation object
   */
  OutputPDB(Node *parent, Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~OutputPDB();

  /*!
   * Open the output stream if we are to write to one file only
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
