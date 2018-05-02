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



#include "gm_scalar_times_scalar_sum.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterScalarTimesScalarSum> grid_meter_scalar_times_scalar_sum("ScalarTimesScalarSum");


/*--- GridMeterScalarTimesScalarSum ---*/

GridMeterScalarTimesScalarSum::GridMeterScalarTimesScalarSum(GridAverager *averager): GridMeter(averager)
{
  init();
}


void GridMeterScalarTimesScalarSum::init()
{
  m_properties.setClassName("ScalarTimesScalarSum");

  m_properties.setDescription("Measure the product of the sum (over particles) of two scalar fields. This Meter simply computes two sums and multiplies them without dividing through the number of particles anywhere.");

  STRINGPC
    (symbol1, m_scalar1_name,
     "Name of the first scalar field to measure.");

  STRINGPC
    (symbol2, m_scalar2_name,
     "Name of the second scalar field to measure.");

//   INTPC
//     (exponent, m_exponent, 0,
//      "Which exponent to measure of the scalar field.");

  m_scalar1_name = "scalar";
  m_scalar2_name = "scalar";

}


void GridMeterScalarTimesScalarSum::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterScalarTimesScalarSum::setup", "Please specify a species.");

  m_result_o = 
    M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + m_scalar1_name + "_time_" + m_scalar2_name + "_sum",
       DataFormat::VECTOR_DOUBLE).offset;

  m_scalar1_o = Particle::s_tag_format[m_colour].attrByName(m_scalar1_name).offset;
  m_scalar2_o = Particle::s_tag_format[m_colour].attrByName(m_scalar2_name).offset;
}


void GridMeterScalarTimesScalarSum::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *results = data->vectorDoubleByOffset(m_result_o).value();
  int n_cells = M_GRID_AVERAGER->nCells();

  m_value1.resize(n_cells);
  m_value2.resize(n_cells);

  results->resize(n_cells);

  for(size_t index = 0; index < n_cells; ++index) {
    m_value1[index] = 0;
    m_value2[index] = 0;
  }

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       
       m_value1[index] += __iSLFE->tag.doubleByOffset(m_scalar1_o); 
       m_value2[index] += __iSLFE->tag.doubleByOffset(m_scalar2_o); 
       
     }
     );

   for(size_t index = 0; index < n_cells; ++index)  
     (*results)[index] = m_value1[index]*m_value2[index];
    
}

