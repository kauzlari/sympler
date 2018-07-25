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



#include "gm_scalar_fourier_cos.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"

GridMeter_Register<GridMeterScalarFourierCos> grid_meter_scalar_fourier_cos("ScalarFourierCos");


/*--- GridMeterScalarFourierCos ---*/

GridMeterScalarFourierCos::GridMeterScalarFourierCos(GridAverager *averager): GridMeterScalar(averager)
{
  init();
}


void GridMeterScalarFourierCos::init()
{
  m_properties.setClassName("ScalarFourierCos");

  m_properties.setDescription("Measure Fourier-cos component of a scalar field.");


  POINTPC
    (kVec, m_kVec,
     "The wave-vector WITHOUT the factor 2Pi which will automatically be appended. ");
      
  for(size_t i = 0; i < SPACE_DIMS; ++i)
    m_kVec[i] = HUGE_VAL;

  //   INTPC(direction, m_dir, -2, "The cartesian direction of the wave described by the wave number given by the attribute 'waveNumber'.");
    
  //   m_dir = -1;

//   DOUBLEPC(waveLength, m_lambda, -HUGE_VAL, "The wave length, of the plane wave vector in the cartesian direction specified by the attribute 'direction'.");

//   m_lambda = HUGE_VAL;

}


void GridMeterScalarFourierCos::setup()
{
  GridMeterScalar::setup();

//   if(m_dir != 0 && m_dir != 1 && m_dir != 2)
//     throw gError("GridMeterScalarFourierCos::setup", "Attribute 'direction' has a forbidden value \"" + ObjToString(m_dir) + "\"! It must have value \"0\", \"1\" or \"2\" representing one of the cartesian directions.");

//   if(m_lambda == HUGE_VAL)
//     throw gError("GridMeterScalarFourierCos::setup","Please define attribute 'waveLength'!");

  bool abort = true;
  for(size_t i = 0; i < SPACE_DIMS; ++i)
    abort = abort && (m_kVec[i] == HUGE_VAL);

  if(abort)
    throw gError("GridMeterScalarFourierCos::setup","k-vector does not seem to be well defined by attribute 'kVec'.");

//   m_k = 2*M_PI/m_lambda;
    
  m_kVec = 2*M_PI*m_kVec;

}


void GridMeterScalarFourierCos::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *moments = data->vectorDoubleByOffset(m_moment_o).value();
  int n_cells = M_GRID_AVERAGER->nCells();

  moments->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       // the argument of the cos is a scalar product
       double e = __iSLFE->tag.doubleByOffset(m_scalar_o)*cos(m_kVec*__iSLFE->r);

       if(m_moment != 1) {
	 double temp = e;
	 for(size_t pow = 2; pow <= m_moment; ++pow)
	   e*=temp;
       }

       /* fixme!!! Slow! Better divide only once per cell. */
       (*moments)[index] += e/n;
     }
    );
}

