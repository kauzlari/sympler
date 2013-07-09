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



#include "gm_scalar_sum.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterScalarSum> grid_meter_scalar_sum("ScalarSum");


/*--- GridMeterScalarSum ---*/

GridMeterScalarSum::GridMeterScalarSum(GridAverager *averager): GridMeter(averager)
{
  init();
}


void GridMeterScalarSum::init()
{
  m_properties.setClassName("ScalarSum");

  m_properties.setDescription("Measure f^exponent of the sum (over particles) of a scalar field f. This Meter simply computes the sum without dividing through the number of particles.");

  STRINGPC
    (symbol, m_scalar_name,
     "Name of the scalar field to measure.");

  INTPC
    (exponent, m_exponent, 0,
     "Which exponent to measure of the scalar field.");

  m_scalar_name = "scalar";
  m_exponent = 1;
}


void GridMeterScalarSum::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterScalarSum::setup", "Please specify a species.");

  m_exponent_o = 
    M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + m_scalar_name + "^" + ObjToString(m_exponent) + "_sum",
       DataFormat::VECTOR_DOUBLE).offset;

  m_scalar_o = Particle::s_tag_format[m_colour].attrByName(m_scalar_name).offset;
}


void GridMeterScalarSum::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *exponents = data->vectorDoubleByOffset(m_exponent_o).value();
  int n_cells = M_GRID_AVERAGER->nCells();

  exponents->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       (*exponents)[index] += __iSLFE->tag.doubleByOffset(m_scalar_o); 
     }
     );

  for(size_t index = 0; index < n_cells; ++index)  
    (*exponents)[index] = pow((*exponents)[index], m_exponent);
    
}

