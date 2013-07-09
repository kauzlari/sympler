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



#include "gm_self_diffusion.h"

#include "phase.h"
#include "controller.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterSelfDiffusion> grid_meter_self_diffusion("SelfDiffusion");


#define M_SIMULATION  ((Simulation*) phase->parent())
#define M_CONTROLLER  M_SIMULATION->controller()


/*--- GridMeterSelfDiffusion ---*/

GridMeterSelfDiffusion::GridMeterSelfDiffusion(GridAverager *averager): GridMeter(averager), m_first_measurement(true)
{
    init();
}


void GridMeterSelfDiffusion::init()
{
  m_properties.setClassName("SelfDiffusion");

  m_properties.setDescription("Measure the self diffusion coefficient. This Meter requires a Symbol for the displacement which is computed elsewhere. ");

  STRINGPC
    (displacement, m_dispName,
     "Name of the displacement to measure.");

  m_dispName = "displacement";

}


void GridMeterSelfDiffusion::setup()
{
  GridMeter::setup();

  m_D_o =
    M_GRID_AVERAGER->format().addAttribute
    (m_species + "_self_diffusion", DataFormat::VECTOR_DOUBLE).offset;

  m_dispOffset = Particle::s_tag_format[m_colour].attrByName(m_dispName).offset;

}


void GridMeterSelfDiffusion::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *sd = data->vectorDoubleByOffset(m_D_o).value();

  double factor = 1/(6*M_CONTROLLER->time());

  //  MSG_DEBUG("GridMeterSelfDiffusion::measure", "...");

  sd->resize(M_GRID_AVERAGER->nCells());

  if (m_first_measurement) {
    m_first_measurement = false;

    /* Reset the diplacement field of each particle */
    /*    MSG_DEBUG
      ("GridMeterSelfDiffusion::measure",
      "Resetting displacement.");*/

    FOR_EACH_FREE_PARTICLE_C
      (phase, m_colour,
       m_initDisp = __iSLFE->tag.pointByOffset(m_dispOffset);
       );
  }

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       /* fixme!!! Slow! Better divide only once per cell. */
       (*sd)[index] += 
	 factor*(__iSLFE->tag.pointByOffset(m_dispOffset) - m_initDisp).absSquare()/n;
     }
    );
}

