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



#include "gm_temperature_ie.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"
#include "integrator_energy.h"

#define M_METER  ((Meter*) m_parent)
#define M_SIMULATION ((Simulation*) M_METER->parent())
#define M_CONTROLLER M_SIMULATION->controller()

GridMeter_Register<GridMeterTemperatureIE> grid_meter_temperature_ie("TemperatureFromInternalEnergy");


/*--- GridMeterTemperatureIE ---*/

GridMeterTemperatureIE::GridMeterTemperatureIE(GridAverager *averager): GridMeter(averager)
{
    init();
}


void GridMeterTemperatureIE::init()
{
  m_properties.setClassName("TemperatureFromInternalEnergy");

  m_properties.setDescription("Measure the local temperature using the internal energy variable.");
}


void GridMeterTemperatureIE::setup()
{
  Integrator *i;

  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterTemperatureIE", "Please specify a species.");

  i = M_CONTROLLER->findIntegrator("IntegratorEnergy", m_species);
  assert(typeid(*i) == typeid(IntegratorEnergy));

  m_ie = (IntegratorEnergy *) i;

  if (!m_ie)
    throw gError("GridMeterTemperatureIE", "I need an IntegratorEnergy for '" + m_species + "'.");

  m_t = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_temperature_mean", DataFormat::VECTOR_DOUBLE).offset;
  m_tv = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_temperature_variance", DataFormat::VECTOR_DOUBLE).offset;

  m_rt = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_reciprocal_temperature_mean", DataFormat::VECTOR_DOUBLE).offset;
  m_rtv = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_reciprocal_temperature_variance", DataFormat::VECTOR_DOUBLE).offset;
}


void GridMeterTemperatureIE::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *temperatures = data->vectorDoubleByOffset(m_t).value();
  vector<double> *tvariances = data->vectorDoubleByOffset(m_tv).value();
  vector<double> *rtemperatures = data->vectorDoubleByOffset(m_rt).value();
  vector<double> *rtvariances = data->vectorDoubleByOffset(m_rtv).value();
  size_t n_cells = M_GRID_AVERAGER->nCells();

  temperatures->resize(n_cells);
  tvariances->resize(n_cells);
  rtemperatures->resize(n_cells);
  rtvariances->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       /* fixme!!! Slow! Better divide only once per cell. */
       double rT = m_ie->reciprocalTemperature(*__iSLFE);
       double T = 1/rT;

       (*temperatures)[index] += T/n;
       (*tvariances)[index] += T*T/n;

       (*rtemperatures)[index] += rT/n;
       (*rtvariances)[index] += rT*rT/n;
     }
    );
}

