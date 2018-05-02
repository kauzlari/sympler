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



#include "gm_density.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterDensity> grid_meter_density("Density");


/*--- GridMeterDensity ---*/

GridMeterDensity::GridMeterDensity(GridAverager *averager)
  : GridMeter(averager)
{
    init();
}


void GridMeterDensity::init()
{
  m_properties.setClassName("Density");

  m_properties.setDescription("Measure the particle density in the grid cells, i.e. particles per cell volume.");
}


void GridMeterDensity::setup()
{
  GridMeter::setup();

  MSG_DEBUG("GridMeterDensity::setup", "m_species = " << m_species);

  if (M_GRID_AVERAGER) {
    m_dm_o = 
      M_GRID_AVERAGER->format().addAttribute
        (m_species + "_density_mean", DataFormat::VECTOR_DOUBLE).offset;
    m_dv_o = 
      M_GRID_AVERAGER->format().addAttribute
        (m_species + "_density_variance", DataFormat::VECTOR_DOUBLE).offset;
  }
}


void GridMeterDensity::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *densities = data->vectorDoubleByOffset(m_dm_o).value();
  vector<double> *variances = data->vectorDoubleByOffset(m_dv_o).value();
  size_t n_cells = M_GRID_AVERAGER->nCells();
  
  densities->resize(n_cells);
  variances->resize(n_cells);

  if (m_colour != ALL_COLOURS) {
    for (size_t i = 0; i < n_cells; i++) {
      double d = M_GRID_AVERAGER->nParticles(m_colour, i)/M_GRID_AVERAGER->volume(i);
      
      (*densities)[i] += d;
      (*variances)[i] += d*d;
//        MSG_DEBUG("GridMeterDensity::measure", "variance" << i << "=" 
//            << (*variances)[i]);    
    }
  } else {
    size_t nc = phase->nColours();

    for (size_t c = 0; c < nc; ++c) {
      for (size_t i = 0; i < n_cells; i++) {
        double d = M_GRID_AVERAGER->nParticles(c, i)/M_GRID_AVERAGER->volume(i);
      
        (*densities)[i] += d;
        (*variances)[i] += d*d;
/*        MSG_DEBUG("GridMeterDensity::measure", "variance" << i << "=" << (*variances)[i]);*/
      }
    }
  }
}

