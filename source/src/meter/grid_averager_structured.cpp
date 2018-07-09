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



#include "grid_averager_structured.h"

#include "simulation.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<GridAveragerStructured> grid_averager_structured("GridAveragerStructured");


#define IDX_TIME 0
#define IDX_N_CELLS 1
#define IDX_CELL_POSITIONS 2
#define IDX_DATA_FORMAT 3


#define INDEX_FROM_POS(i, j, k)  (k*m_n.y+j)*m_n.x+i


#define M_BOUNDARY   ((Simulation*) m_parent)->phase()->boundary()


//---- Constructors/Destructor ----

GridAveragerStructured::GridAveragerStructured(Simulation *simulation)
  : GridAverager(simulation)/*, m_step(0)*/
{
  init();
}


GridAveragerStructured::~GridAveragerStructured()
{
}



//---- Methods ----

void GridAveragerStructured::init()
{
  /* Properties */
  m_properties.setClassName("GridAveragerStructured");

  m_properties.setDescription(
    "Determines local properties on a structured 3D grid. "
    "The properties "
    "may be added as modules. Read the help for 'GridMeters' for more information."
    "The output must be produced by Postprocessors which "
    "may be added as modules. Read the help for 'Postprocessors' for more information."
  );

  for (int i = 0; i < SPACE_DIMS; i++) {
    m_properties.addProperty
      ("n" + string(1, 'X'+i), PropertyList::INT, &m_n[i],
       new PLCIntGreater(0),
       "Number of grid cells in " + string(1, 'x'+i) + "-"
       "direction.");
    m_n[i] = 10;
  }    

  /* Format for the output pipe */
  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute(IDDF_N_CELLS, DataFormat::INT_POINT);
  m_format.addAttribute("cell_positions", DataFormat::VECTOR_POINT);
  m_format.addAttribute(IDDF_DATA_FORMAT, DataFormat::INT);

  m_idx_data_start = IDX_DATA_FORMAT+1;
}


void GridAveragerStructured::setup()
{
  m_n_cells = m_n.x*m_n.y*m_n.z;

  GridAverager::setup();
}


void GridAveragerStructured::aboutToStart()
{
  m_volume = ((Simulation*) m_parent)->phase()->boundary()->cuboidVolume()/m_n_cells;

  MSG_DEBUG("GridAveragerStructured::setup", "cell volume = " << m_volume);

  /* Make sure, setup of the postprocessor is called *last*. */
  GridAverager::aboutToStart();
}


/*virtual*/ void GridAveragerStructured::measureNow(const double& time)
{
  if (!m_step) {
    m_data = m_format.newData();

    m_data->doubleByIndex(IDX_TIME) = time;
    m_data->intPointByIndex(IDX_N_CELLS) = m_n;
    m_data->vectorPointByIndex(IDX_CELL_POSITIONS)->resize(m_n_cells);
    m_data->intByIndex(IDX_DATA_FORMAT) = DF_STRUCTURED_GRID;

    vector<point_t> *positions = m_data->vectorPointByIndex(IDX_CELL_POSITIONS).value();
    point_t offset, size;

    offset = M_BOUNDARY->boundingBox().corner1;
    size = M_BOUNDARY->boundingBox().size();
       
    for (int k = 0; k < m_n.z; ++k) {
      for (int j = 0; j < m_n.y; ++j) {
        for (int i = 0; i < m_n.x; ++i) {
          point_t p;

          p.x = (i+0.5)*size.x/m_n.x;
          p.y = (j+0.5)*size.y/m_n.y;
          p.z = (k+0.5)*size.z/m_n.z;

          (*positions)[INDEX_FROM_POS(i, j, k)] = p+offset;
/*         MSG_DEBUG("GridAveragerStructured::measureNow", "working on cell with index " 
           << INDEX_FROM_POS(i, j, k) << " and pos=" << p+offset);*/
        }
      }
    }
  }

  findLocations();
  average(m_data);
  ++m_step;

  if (m_step == m_avg_steps) {
    finishStep(m_data);
    m_step = 0;
  }
}


void GridAveragerStructured::findLocations()
{
  Phase *phase = ((Simulation*) m_parent)->phase();
  point_t offset, size;

  offset = phase->boundary()->boundingBox().corner1;
  size = phase->boundary()->boundingBox().size();

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
     int_point_t cell;

     for (int j = 0; j < SPACE_DIMS; j++) {
       cell[j] = (int) ((__iSLFE->r[j]-offset[j])*m_n[j]/size[j]);
     }
     int index = INDEX_FROM_POS(cell.x, cell.y, cell.z);
     
     if(index < 0 || index > m_n_cells) 
       throw gError("GridAveragerStructured::findLocations", "Invalid index " + ObjToString(index) + " for particle " + ObjToString(__iSLFE->mySlot) + " with coordinates (" + ObjToString(__iSLFE->r.x) + ", " + ObjToString(__iSLFE->r.y) + ", " + ObjToString(__iSLFE->r.z) + ")");
     
     m_locations[c][__iSLFE->mySlot] = index;
     ++(m_n_particles[c][index]);
     ++(m_n_total_particles[index]);
    );
}


void GridAveragerStructured::flush()
{
  if (m_step)
    finishStep(m_data);

  GridAverager::flush();
}

