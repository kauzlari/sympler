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



#include "gm_temperature.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterTemperature> grid_meter_temperature("Temperature");


/*--- GridMeterTemperature ---*/

GridMeterTemperature::GridMeterTemperature(GridAverager *averager): GridMeter(averager)
{
    init();
}


void GridMeterTemperature::init()
{
  m_properties.setClassName("Temperature");

  m_properties.setDescription("Measure the local temperature.");
}


void GridMeterTemperature::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterTemperature", "Please specify a species. Measurements over all species are currently not supported by this GridMeter.");

  m_TIndex = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_temperature", DataFormat::VECTOR_DOUBLE).index;
  m_vCMIndex = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_vCM", DataFormat::VECTOR_POINT).index;
}


void GridMeterTemperature::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *temperatures = data->vectorDoubleByIndex(m_TIndex).value();
  vector<point_t> *vCM = data->vectorPointByIndex(m_vCMIndex).value();
  size_t n_cells = M_GRID_AVERAGER->nCells();

  temperatures->resize(n_cells);
  vCM->resize(n_cells);
  m_temperatures.resize(n_cells);
  m_vCM.resize(n_cells);

  for(size_t cell = 0; cell < n_cells; ++cell)
  {
    m_temperatures[cell] = 0;
    for (int i = 0; i < SPACE_DIMS; i++)
      ((m_vCM)[cell][i]) = 0;
	
  }
  // 3/2 kT = 1/2 m v^2  => T = v^2/3 with k=m=1
  // we do <(v-<v>)^2> = <v^2>-<v>^2

  // first we sum up in the loop over particles
  FOR_EACH_FREE_PARTICLE_C
    (phase,
     m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);

     if (M_GRID_AVERAGER->nParticles(c, index)) {
       m_temperatures[index] += __iSLFE->v.absSquare();
       m_vCM[index] += __iSLFE->v;
     }
    );

  size_t nOfCols = phase->nColours();
  for(size_t cell = 0; cell < n_cells; ++cell)
  {

    size_t N = M_GRID_AVERAGER->nParticles(m_colour, cell);

    if(N)
      {
	// this is now <v^2>/3
	(m_temperatures)[cell] /= 3*N;
	// here we do <v> 
	for (int i = 0; i < SPACE_DIMS; i++)
	  (m_vCM)[cell][i] /= N;
	// now the difference <v^2>-<v>^2, where <v>^2/3 is done on the fly
	(m_temperatures)[cell] -= (m_vCM)[cell].absSquare()/3;
	// transfer the results to the grid_averager
	(*temperatures)[cell] += (m_temperatures)[cell];
	(*vCM)[cell] += m_vCM[cell]; 
      }
  }
}

