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



#ifndef __OUTPUT_FILE_H
#define __OUTPUT_FILE_H

#include <list>
#include <fstream>

#include "output.h"

using namespace std;


/*!
 * Writes data to a text file
 */
class OutputFile: public Output
{
protected:
  /*!
   * Format description of the data passed on by the parent
   * \a Postprocessor
   */
  DataFormat *m_input_format;

  /*!
   * The columns from the parent data to be written
   */
  list<int> m_columns;

  /*!
   * The filename to write to
   */
  string m_filename;

  /*!
   * String identifier for the columns to be written. Columns are 
   * seperated by the pipe '|' symbol.
   */
  string m_cols;

  /*!
   * Start of a commentary block (for the header)
   */
  string m_comment_start;

  /*!
   * End of a commentary block (for the header)
   */
  string m_comment_end;

  /*!
   * The output stream
   */
  ofstream m_s;

  /*!
   * Write a header (encapsulated in a commentary block)?
   */
  bool m_do_head;

  /*!
   * Write each data that is pushed to this output module into
   * a seperate file?
   */
  bool m_multiple_files; 

  /*!
   * If set to false, ascii files are produced, if set to true, binary files. In binary mode, of course no header is written (\a m_do_head has no effect).
   */
  bool m_binary;

  /*!
   * Are all columns vectors (i.e. VECTOR_DOUBLE, VECTOR_POINT, ...)?
   * If this is the case, and multiple files should be written,
   * the columns are written as columns into the file (because they
   * can be nicely aligned).
   * Irrelevant if \a m_binary == true
   */
  bool m_all_vectors;

  /*!
   * The current output step (will be appended to the file name 
   * if multiple files are to be written).
   */
  int m_step_counter;
    
  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Write the file header
   */
  virtual void writeHeader();

  /*!
   * Open the stream for writing
   */
  virtual void openStream();

  /*!
   * Close the stream
   */
  virtual void closeStream();
    
public:
  /*!
   * Constructor
   * @param parent Parent node in the \a Postprocessor chain
   * @param simulation Pointer to the main \a Simulation object
   */
  OutputFile(Node *parent, Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~OutputFile();

  /*!
   * Parse the columns string \a m_cols and open file for writing
   * if we are writing to just one file
   */
  virtual void setup();

  /*!
   * Additional setup
   */
  virtual void setupAfterParticleCreation();
  
  /*!
   * Called by the parent to pass on its \a DataFormat
   */
  virtual void describeInput(DataFormat *input_format);

  /*!
   * Data from the parent \a Postprocessor or \a Meter
   */
  virtual void push(data_sp data);

  /*!
   * Close the file if we're just writing to one file
   */
  virtual void flush();
};

#endif
