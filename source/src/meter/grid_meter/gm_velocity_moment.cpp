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



#include "gm_velocity_moment.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterVelocityMoment> grid_meter_velocity_moment("VelocityMoment");


/*--- GridMeterVelocityMoment ---*/

GridMeterVelocityMoment::GridMeterVelocityMoment(GridAverager *averager)
  : GridMeter(averager)
{
  init();
}


void GridMeterVelocityMoment::init()
{
  m_properties.setClassName("VelocityMoment");

  m_properties.setDescription
    ("Maximum likelyhood estimation for the moments of the velocity distribution.");

  INTPC
    (order, m_order, 0,
     "Order of the moment to determine, i.e. <v> -> order = 1, <v^2> -> order = 2 ...");
}


void GridMeterVelocityMoment::setup()
{
  GridMeter::setup();

  if (M_GRID_AVERAGER) {
    m_offset = M_GRID_AVERAGER->format().addAttribute
      (m_species + "_vavg^" + ObjToString(m_order), DataFormat::VECTOR_DOUBLE).offset;

    m_offset_by_dir = M_GRID_AVERAGER->format().addAttribute
      (m_species + "_vxyz^" + 
       ObjToString(m_order), DataFormat::VECTOR_POINT).offset;
  }
}


void GridMeterVelocityMoment::measure(data_sp data)
{
  point_t null_point = {{{ 0, 0, 0 }}};

  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();

  vector<double> *avg_velocities = data->vectorDoubleByOffset(m_offset).value();
  vector<point_t> *velocities = data->vectorPointByOffset(m_offset_by_dir).value();

  size_t n_cells = M_GRID_AVERAGER->nCells();

  avg_velocities->clear();
  avg_velocities->resize(n_cells, 0);

  velocities->clear();
  velocities->resize(n_cells, null_point);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     point_t p;
       
     int cell = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, cell);

     for (int j = 0; j < SPACE_DIMS; ++j)
       p[j] = pow(__iSLFE->v[j], m_order);
       
     (*avg_velocities)[cell] += (p.x+p.y+p.z)/(3*n);
     (*velocities)[cell] += p/n;
    );
}

