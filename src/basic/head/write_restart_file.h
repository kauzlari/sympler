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



#ifndef __WRITE_RESTART_FILE_H
#define __WRITE_RESTART_FILE_H 

#include "callable.h"
#include "integrator_energy.h"

using namespace std;

/* --- WriteRestartFile --- */

/*!
 * Write a restart file during the simulation
 */
class WriteRestartFile : public Callable
{
protected:
  /*!
   * File name of the restart file
   */
  string m_file_name;
 
  /*!
   * Interval in which to write the restart file
   */
  size_t m_write_every;

  /*!
   * Current step
   */
  size_t m_step;

  /*!
   * Current file number to write
   */
  size_t m_file_no;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  WriteRestartFile(Simulation* simulation);

  /*!
   * Destructor
   */
  virtual ~WriteRestartFile() {}

  /*!
   * Write the restart file if \a m_write_every steps have passed.
   */
  virtual void call(size_t timestep);

  /*!
   * Initialize \a m_step and \a m_file_no
   */
  virtual void setup();
};

#endif
