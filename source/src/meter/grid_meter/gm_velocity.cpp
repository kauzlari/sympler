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



#include "gm_velocity.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterVelocity> grid_meter_velocity("Velocity");


/*--- GridMeterVelocity ---*/

GridMeterVelocity::GridMeterVelocity(GridAverager *averager): GridMeter(averager)
{
    init();
}


void GridMeterVelocity::init()
{
  m_properties.setClassName("Velocity");

  m_properties.setDescription("Measure the average velocity vector.");
}


void GridMeterVelocity::setup()
{
  GridMeter::setup();
  MSG_DEBUG("GridMeterVelocity::setup", "m_species=" << m_species);
  m_vm_o =
    M_GRID_AVERAGER->format().addAttribute
      (m_species + "_velocity_mean", DataFormat::VECTOR_POINT).offset;
  m_vv_o =
    M_GRID_AVERAGER->format().addAttribute
      (m_species + "_velocity_variance", DataFormat::VECTOR_POINT).offset;

  m_avm_o =
    M_GRID_AVERAGER->format().addAttribute
      (m_species + "_abs_velocity_mean", DataFormat::VECTOR_DOUBLE).offset;
  m_avv_o =
    M_GRID_AVERAGER->format().addAttribute
      (m_species + "_abs_velocity_variance", DataFormat::VECTOR_DOUBLE).offset;
}


void GridMeterVelocity::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<point_t> *velocities = data->vectorPointByOffset(m_vm_o).value();
  vector<point_t> *variances = data->vectorPointByOffset(m_vv_o).value();    
  vector<double> *avelocities = data->vectorDoubleByOffset(m_avm_o).value();
  vector<double> *avariances = data->vectorDoubleByOffset(m_avv_o).value();

  velocities->resize(M_GRID_AVERAGER->nCells());
  variances->resize(M_GRID_AVERAGER->nCells());
  avelocities->resize(M_GRID_AVERAGER->nCells());
  avariances->resize(M_GRID_AVERAGER->nCells());

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       /* fixme!!! Slow! Better divide only once per cell. */
       (*velocities)[index] += __iSLFE->v/n;
       (*avelocities)[index] += __iSLFE->v.abs()/n;
       (*avariances)[index] += __iSLFE->v.absSquare()/n;

       for (int j = 0; j < SPACE_DIMS; j++)
         (*variances)[index][j] += __iSLFE->v[j]*__iSLFE->v[j]/n;
     }
    );
}

