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



#include "grid_meter.h"

#include "simulation.h"
#include "manager_cell.h"
#include "grid_averager.h"

REGISTER_SMART_ENUM
(GridMeter_Factory,
 "GridMeters measure physical quantities in grid cells. If not stated differently or not done per definition, the measured quantitiy is measured per particle in the grid cell."
);

#define M_AVERAGER  ((GridAverager*) m_parent)
#define M_SIMULATION ((Simulation*) M_AVERAGER->parent())
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


/*--- GridMeter ---*/

GridMeter::GridMeter(GridAverager *averager): Node(averager), m_colour(ALL_COLOURS)
{
  init();
}


GridMeter::~GridMeter()
{
}


void GridMeter::init()
{
  STRINGPC
    (species, m_species,
     "Defines the species for the measurement. Some"
     " GridMeters can measure either all or exactly one species.");

  m_species = "UNDEF";
}


void GridMeter::setup()
{
  Node::setup();

  if (m_species != "UNDEF")
    m_colour = M_MANAGER->getColour(m_species);
}

/* --- Utitlity functions --- */

void GridMeter::calcVariance(const vector_double_sp &values, vector_double_sp &values_sq, size_t n_steps)
{
  assert(values->size() == values_sq->size());

  int n_cells = values->size();

  // we have to divide once through the number of timesteps because the RHS is 
  // proportional to the square of, so to divisions are needed, and the GridAverager 
  // makes only one
  for (int i = 0; i < n_cells; i++) {
    (*values_sq)[i] -= (*values)[i]*(*values)[i]/n_steps;
  }
}


void GridMeter::calcVariance(const vector_point_sp &values, vector_point_sp &values_sq, size_t n_steps)
{
  assert(values->size() == values_sq->size());

  int n_cells = values->size();

  for (int i = 0; i < n_cells; i++) {
    for (int j = 0; j < SPACE_DIMS; j++) {
      (*values_sq)[i][j] -= (*values)[i][j]*(*values)[i][j]/n_steps;
    }
  }
}


