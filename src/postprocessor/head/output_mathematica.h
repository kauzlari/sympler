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



#ifndef __OUTPUT_MATHEMATICA_H
#define __OUTPUT_MATHEMATICA_H

#include "output_file.h"

/*!
 * Output file in Mathematica notation (i.e., use curly brackets)
 */
class OutputMathematica: public OutputFile
{
 protected:
  /*!
   * If the output is written to one file, we need to encapsulate
   * the data in an extra "{ ... }". This variable tells the output to put a
   * comma at the end of each data pushed to this output module.
   */
  bool m_comma;

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
  OutputMathematica(Node *parent, Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~OutputMathematica();

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
