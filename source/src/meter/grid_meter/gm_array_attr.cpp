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



#include "gm_array_attr.h"

#ifdef WITH_ARRAY_TYPES


#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterArrayAttr> grid_meter_array_attr("Array");


/*--- GridMeterScalar ---*/

GridMeterArrayAttr::GridMeterArrayAttr(GridAverager *averager): GridMeter(averager)
{
  init();
}


void GridMeterArrayAttr::init()
{
  m_properties.setClassName("Array");

  m_properties.setDescription("Output a vector or matrix of MArray2D format.");

  STRINGPC
      (symbol, m_name,
       "Name of the vector or matrix.");

  m_name = "stateVector";
}


void GridMeterArrayAttr::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterArrayAttr::setup", "Please specify a species.");

  m_output_o = 
      M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + m_name,
       DataFormat::MArray2D).offset;

  m_offset = Particle::s_tag_format[m_colour].attrByName(m_name).offset;
  MSG_DEBUG("GridMeterArrayAttr::setup", "m_offset = " << m_offset << ", m_name = " << m_name);
}


void GridMeterArrayAttr::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  
  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
     
     // FIXME!: the main authors of this file (cenovai, lieneman) should check if the next lines 
     // make sense and work; it seems to compile at least
     // there is no operator+= at the moment
     (*(data->array2dDoubleByOffset(m_output_o))) = (*(data->array2dDoubleByOffset(m_output_o))) + (__iSLFE->tag.array2dDoubleByOffset(m_offset)).deepCopy()->scalarmult(1./n);
     // 			   data->array2dDoubleByOffset(m_output_o) += (__iSLFE->tag.array2dDoubleByOffset(m_offset)).deepCopy()->scalarmult(1./n);
     );
}

#endif
