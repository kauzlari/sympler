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



#include "wf_tabulated.h"

#include "simulation.h"

const WeightingFunction_Register<TabulatedWeightingFunctionWithWall> 
  tabulated_weighting_function_with_wall("TabulatedWeightingFunctionWithWall");


#define M_SIMULATION ((Simulation*) m_parent)


/* --- TabulatedWeightingFunctionWithWall --- */

/* Constructor/Destructor */


TabulatedWeightingFunctionWithWall::TabulatedWeightingFunctionWithWall(Node *parent): WeightingFunctionWithWall(parent)
{
  init();
}

TabulatedWeightingFunctionWithWall::~TabulatedWeightingFunctionWithWall()
{
}


/* Methods */

void TabulatedWeightingFunctionWithWall::setup() {
  point_t g = {{{ 0, 0, 0 }}};

  MSG_DEBUG("TabulatedWeightingFunctionWithWall::setup", "Calculation weighting tables.");

  m_wf = M_SIMULATION->findWeightingFunction(m_wf_name);

  cout << m_wf << endl;

  if (!m_wf)
    throw gError
      ("TabulatedWeightingFunctionWithWall::setup",
       "Please specify a weighting function to discretize.");

  m_cutoff = m_wf->cutoff();

  cout << "0.1" << endl;

  m_interpolate_table.setSize(m_n_bins+1, m_n_bins);
  m_weight_table.setSize(m_n_bins+1, m_n_bins);
  m_local_gradient_table.setSize(m_n_bins+1, m_n_bins);

  cout << "1" << endl;

  for (size_t j = 0; j < m_n_bins; j++) {
    for (size_t i = 0; i < m_n_bins; i++) {
      m_interpolate_table(i, j) = 
	m_wf->interpolateWithDist(i*m_cutoff/m_n_bins, ((double) j)/(m_n_bins-1));
      m_weight_table(i, j) =
	m_wf->weightWithDist(i*m_cutoff/m_n_bins, ((double) j)/(m_n_bins-1));
      m_local_gradient_table(i, j) =
	m_wf->localGradientWithDist(i*m_cutoff/m_n_bins, ((double) j)/(m_n_bins-1));
    }

    m_interpolate_table(m_n_bins, j) = 0;
    m_weight_table(m_n_bins, j) = 0;
    m_local_gradient_table(m_n_bins, j) = g;
  }
}

void TabulatedWeightingFunctionWithWall::init() {
  m_properties.setClassName("TabulatedWeightingFunctionWithWall");
  m_properties.setDescription(
      "Discretises a given weighting function into a lookup table. Currently, the discretisation is evenly spaced.");

  STRINGPC
    (weightingFunction, m_wf_name,
     "Name of the weighting function to be discretised.");

  INTPC
    (nBins, m_n_bins, 0,
     "Number of bins for the discretisation of the weighting functions.");

  m_wf_name = "default";
  m_n_bins = 100;
}

