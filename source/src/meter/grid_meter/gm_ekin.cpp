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



#include "gm_ekin.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterEKin> grid_meter_ekin("KineticEnergy");


/*--- GridMeterEKin ---*/

GridMeterEKin::GridMeterEKin(GridAverager *averager): GridMeter(averager)
{
  init();
}


void GridMeterEKin::init()
{
  m_properties.setClassName("KineticEnergy");

  m_properties.setDescription("Measure the kinetic energy (assuming a particle mass of 1).");
}


void GridMeterEKin::setup()
{
  GridMeter::setup();

  if (M_GRID_AVERAGER)
    m_offset = M_GRID_AVERAGER->format().addAttribute
      (m_species + "_kinetic_energy", DataFormat::VECTOR_DOUBLE).offset;
}


void GridMeterEKin::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *energies = data->vectorDoubleByOffset(m_offset).value();
  size_t n_cells = M_GRID_AVERAGER->nCells();

  energies->clear();
  energies->resize(n_cells, 0);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
       
     int cell = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, cell);
       
     if (n) {
       (*energies)[cell] += __iSLFE->v.absSquare()/(2*n);
     }
    );
}

