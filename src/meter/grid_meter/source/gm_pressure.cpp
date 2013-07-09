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



#include "gm_pressure.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterPressure> grid_meter_pressure("Pressure");

#define M_METER  ((Meter*) m_parent)
#define M_SIMULATION ((Simulation*) M_METER->parent())
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()


/*--- GridMeterPressure ---*/

GridMeterPressure::GridMeterPressure(GridAverager *averager): GridMeter(averager)
{
    init();
}


void GridMeterPressure::init()
{
  m_properties.setClassName("Pressure");

  m_properties.setDescription
    ("Measure the virial part of the local pressure in the grid cell. By definition this is not a per-particle result.");
    
  STRINGPC
  (stress, m_stress_name,
   "Name of the stress field.");
   
  m_stress_name = "stress";
}


void GridMeterPressure::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterPressure", "Please specify a species.");

  m_p_offset = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_pressure_mean", DataFormat::VECTOR_DOUBLE).offset;
  //  m_pv = M_GRID_AVERAGER->format().addAttribute
  //    (m_species + ":pressure:variance", DataFormat::VECTOR_DOUBLE).offset;
  
  
  m_stress_offset = Particle::s_tag_format[m_colour].attrByName(m_stress_name).offset;
}


void GridMeterPressure::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *p = data->vectorDoubleByOffset(m_p_offset).value();
  //  vector<double> *pv = data->vectorDoubleByOffset(m_pv).value();
  size_t n_cells = M_GRID_AVERAGER->nCells();
//   double factor = 1/(SPACE_DIMS * M_PHASE->cuboidVolume());

  p->resize(n_cells);
  //  pv->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
     
     if (n) {
       double pr = 0.5*(__iSLFE->tag.tensorByOffset(m_stress_offset)).trace()/
             (SPACE_DIMS * M_GRID_AVERAGER->volume(index));
       (*p)[index] += pr;
     }
    );
}

