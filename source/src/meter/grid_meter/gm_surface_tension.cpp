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



#include "gm_surface_tension.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterSurfaceTension> grid_meter_surface_tension("SurfaceTension");

#define M_METER  ((Meter*) m_parent)
#define M_SIMULATION ((Simulation*) M_METER->parent())
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()


/*--- GridMeterSurfaceTension ---*/

GridMeterSurfaceTension::GridMeterSurfaceTension(GridAverager *averager): GridMeter(averager)
{
    init();
}


void GridMeterSurfaceTension::init()
{
  m_properties.setClassName("SurfaceTension");

  m_properties.setDescription
    ("Measure the surface tension of an interface. Note: The value returned has "
     "to be multiplied by the width of the slice that is regarded for the pressure "
     "calculation. And the value is not a per-particle value.");

  INTPC
    (perpDir, m_perp_dir, -1,
     "Direction perpendicular to the interface: 0 = x, 1 = y, 2 = z.");
     
  STRINGPC
  (stress, m_stress_name,
   "Name of the stress field.");
   
  m_stress_name = "stress";
}


void GridMeterSurfaceTension::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterSurfaceTension", "Please specify a species.");

  m_st_offset = M_GRID_AVERAGER->format().addAttribute
    (m_species + "_surface_tension", DataFormat::VECTOR_DOUBLE).offset;

  m_perp_dir %= 3;
  m_other_dir1 = (m_perp_dir+1)%3;
  m_other_dir2 = (m_perp_dir+2)%3;
  
  m_stress_offset = Particle::s_tag_format[m_colour].attrByName(m_stress_name).offset;
}


void GridMeterSurfaceTension::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *st = data->vectorDoubleByOffset(m_st_offset).value();
  size_t n_cells = M_GRID_AVERAGER->nCells();

  st->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
     double volume = M_GRID_AVERAGER->volume(index);

     if (n) {
       double Pxx;
       double Pyy;
       double Pzz;

       tensor_t* tensor_stress = &__iSLFE->tag.tensorByOffset(m_stress_offset);  

       Pxx = __iSLFE->v[m_other_dir1]*__iSLFE->v[m_other_dir1] + 0.5*(*tensor_stress)(m_other_dir1, m_other_dir1);
       Pyy = __iSLFE->v[m_other_dir2]*__iSLFE->v[m_other_dir2] + 0.5*(*tensor_stress)(m_other_dir2, m_other_dir2);
       Pzz = __iSLFE->v[m_perp_dir]*__iSLFE->v[m_perp_dir] + 0.5*(*tensor_stress)(m_perp_dir, m_perp_dir);

       (*st)[index] += (Pzz-0.5*(Pxx+Pyy))/volume;
     }
     );
}

