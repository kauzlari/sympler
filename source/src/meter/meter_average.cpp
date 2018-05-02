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



#include "meter_average.h"

#include "simulation.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterAverage> meter_average("MeterAverage");


#define IDX_TIME 0



//---- Constructors/Destructor ----

MeterAverage::MeterAverage(Simulation *simulation)
    : GridAverager(simulation)
{
    init();
}


MeterAverage::~MeterAverage()
{
}



//---- Methods ----

void MeterAverage::init()
{
  /* Properties */
  m_properties.setClassName("MeterAverage");

  m_properties.setDescription(
    "Determines averages over the whole simulation box."
    " For its available modules see the help for 'GridMeters'"
  );

  /* Format for the output pipe */
  m_format.addAttribute("time", DataFormat::DOUBLE);

  m_idx_data_start = IDX_TIME+1;
  m_n_cells = 1;
}


void MeterAverage::setup()
{
  m_volume = ((Simulation*) m_parent)->phase()->boundary()->cuboidVolume();
  m_n_cells = 1;

  GridAverager::setup();
}


/*virtual*/ void MeterAverage::measureNow(const double& time)
{
  if (!m_step) {
    m_data = m_format.newData();

    m_data->doubleByIndex(IDX_TIME) = time;
  }

  findLocations();
  average(m_data);
  m_step++;

  if (m_step == m_avg_steps) {
    finishStep(m_data);
    m_step = 0;
  }
}


void MeterAverage::findLocations()
{
  Phase *phase = ((Simulation*) m_parent)->phase();

  size_t nOfC = phase->nColours();
  for(size_t colour = 0; colour < nOfC; ++colour)
  {
    m_n_particles[colour].clear();
    m_n_particles[colour].push_back(phase->returnNofPartC(colour));
    /* We have only one cell, the whole box */
    if (m_locations[colour].capacity() != phase->returnArrayCapacity(colour))
      m_locations[colour].resize(phase->returnArrayCapacity(colour), 0);
  }    


}


void MeterAverage::flush()
{
    if (m_step)
        finishStep(m_data);

    GridAverager::flush();
}

