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



#include "grid_averager_circular.h"

#include "simulation.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<GridAveragerCircular> grid_averager_circular("GridAveragerCircular");


#define IDX_TIME 0
#define IDX_N_CELLS 1
#define IDX_CELL_POSITIONS 2
#define IDX_DATA_FORMAT 3


#define M_BOUNDARY   ((Simulation*) m_parent)->phase()->boundary()


//---- Constructors/Destructor ----

GridAveragerCircular::GridAveragerCircular(Simulation *simulation)
  : GridAverager(simulation)
{
  init();
}


GridAveragerCircular::~GridAveragerCircular()
{
}



//---- Methods ----

void GridAveragerCircular::init()
{
  /* Properties */
  m_properties.setClassName("GridAveragerCircular");

  m_properties.setDescription(
    "Determines local properties on a structured 3D grid. The properties may be added as modules. "
    "Read the help for 'GridMeters' for more information"
  );

  DOUBLEPC
    (innerRadius, m_inner_radius, 0,
     "Inner radius of the cylindrical geometry to average over.");

  DOUBLEPC
    (outerRadius, m_outer_radius, 0,
     "Outer radius of the cylindrical geometry to average over.");

  INTPC
    (nZ, m_n_z, 0,
     "Number of cells in z-direction. As for now, only one cell is possible.");

  INTPC
    (nR, m_n_r, 0,
     "Number of cells in r-direction.");

  INTPC
    (nPhi, m_n_phi, 0,
     "Number of cells in phi-direction. As for now, only one cell is possible.");


  m_inner_radius = 1;
  m_outer_radius = 10;
  m_n_z = 1;
  m_n_r = 10;
  m_n_phi = 1;


  /* Format for the output pipe */
  m_format.addAttribute("time", DataFormat::DOUBLE, false);
//  m_format.addAttribute(IDDF_N_CELLS, DataFormat::INT_POINT);
//  m_format.addAttribute("cell_positions", DataFormat::VECTOR_POINT);
//  m_format.addAttribute(IDDF_DATA_FORMAT, DataFormat::INT);

//  m_idx_data_start = IDX_DATA_FORMAT+1;
  m_idx_data_start = 1;
}


void GridAveragerCircular::setup()
{
  if (m_n_z != 1)
    throw gError
      ("GridAveragerCircular::setup",
       "As for now, it is only possible to define one cell in z-direction.");
  if (m_n_phi != 1)
    throw gError
      ("GridAveragerCircular::setup",
       "As for now, it is only possible to define one cell in phi-direction.");

  m_n_cells = m_n_z*m_n_r*m_n_phi;

  GridAverager::setup();
}


void GridAveragerCircular::aboutToStart()
{
  point_t size = ((Simulation*) m_parent)->phase()->boundary()->boundingBox().size();
  
  m_diff_r = (m_outer_radius-m_inner_radius)/m_n_cells;

  m_volume.resize(m_n_cells);

  for (int i = 0; i < m_n_cells; i++) {
    m_volume[i] = (M_PI*pow(m_inner_radius + (i+1)*m_diff_r, 2) - M_PI*pow(m_inner_radius + i*m_diff_r, 2))*size.z;
  }

  /* Make sure, setup of the postprocessor is called *last*. */
  GridAverager::aboutToStart();
}


/*virtual*/ void GridAveragerCircular::measureNow(const double& time)
{
  if (!m_step) {
    m_data = m_format.newData();

    m_data->doubleByIndex(IDX_TIME) = time;
//    m_data->intPointByIndex(IDX_N_CELLS) = m_n;
//    m_data->vectorPointByIndex(IDX_CELL_POSITIONS)->resize(m_n_cells);
//    m_data->intByIndex(IDX_DATA_FORMAT) = DF_STRUCTURED_GRID;

/*
    vector<point_t> *positions = m_data->vectorPointByIndex(IDX_CELL_POSITIONS).value();
    point_t offset, size;

    offset = M_BOUNDARY->boundingBox().corner1;
    size = M_BOUNDARY->boundingBox().size();
       
    for (int k = 0; k < m_n.z; k++) {
      for (int j = 0; j < m_n.y; j++) {
        for (int i = 0; i < m_n.x; i++) {
          point_t p;

          p.x = (i+0.5)*size.x/m_n.x;
          p.y = (j+0.5)*size.y/m_n.y;
          p.z = (k+0.5)*size.z/m_n.z;

          (*positions)[INDEX_FROM_POS(i, j, k)] = p+offset;
        }
      }
    }
*/
  }

  findLocations();
  average(m_data);
  m_step++;

  if (m_step == m_avg_steps) {
    finishStep(m_data);
    m_step = 0;
  }
}


void GridAveragerCircular::findLocations()
{
  Phase *phase = ((Simulation*) m_parent)->phase();
//   point_t offset, size;

  size_t nOfC = m_n_particles.size();
  assert(nOfC == phase->nColours());

  for(size_t c = 0; c < nOfC ; ++c) {
    m_n_particles[c].clear();
    m_n_particles[c].resize(m_n_cells, 0);

    if (m_locations[c].capacity() != phase->returnArrayCapacity(c))
      m_locations[c].resize(phase->returnArrayCapacity(c));
  }

  m_n_total_particles.resize(m_n_cells, 0);

  FOR_EACH_FREE_PARTICLE
    (phase,
     int index = (int) ((__iSLFE->r.abs() - m_inner_radius)/m_diff_r);

     if (index < 0)
       index = 0;
     else if (index >= m_n_cells)
       index = m_n_cells-1;
         
     m_locations[c][__iSLFE->mySlot] = index;
     m_n_particles[c][index]++;
     m_n_total_particles[index]++;
    );
}


void GridAveragerCircular::flush()
{
  if (m_step)
    finishStep(m_data);

  GridAverager::flush();
}

