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



#include "write_restart_file.h"


#include "phase.h"
#include "simulation.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()


const Callable_Register<WriteRestartFile> write_restart_file("WriteRestartFile");


WriteRestartFile::WriteRestartFile(Simulation* sim)
  : Callable(sim)
{
  init();
}


void WriteRestartFile::init()
{
  m_properties.setClassName("WriteRestartFile");

  m_properties.setDescription(
    "Output a restart file more often than just at the end of the simulation."
  );

  STRINGPCOUF
    (name, m_file_name,
     "Name of the restart file.");

  INTPC
    (writeEvery, m_write_every, 0,
     "Interval between writing restart files.");

  m_file_name = "restart.pos";
  m_write_every = 100;
}


void WriteRestartFile::setup()
{
  Callable::setup();

  m_step = 0;
  m_file_no = 0;
}


void WriteRestartFile::call(size_t timestep)
{
  Phase *phase = M_PHASE;

  m_step++;

  if (m_step == m_write_every) {
    MSG_DEBUG("WriteRestartFile::call", "Outputting restart file.");

    phase->writeRestartFile(make_filename(m_file_name, m_file_no));
    
    m_file_no++;
    m_step = 0;
  }
}


