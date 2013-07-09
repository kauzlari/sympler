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



#include "gm_local_density.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"
#include "val_calculator_rho.h"

GridMeter_Register<GridMeterLocalDensity> grid_meter_local_density("LocalDensity");

#define STR_LOCAL_DENSITY "local_density"

#define M_AVERAGER  ((GridAverager*) m_parent)
#define M_SIMULATION ((Simulation*) M_AVERAGER->parent())
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


/*--- GridMeterLocalDensity ---*/

GridMeterLocalDensity::GridMeterLocalDensity(GridAverager *averager): GridMeter(averager)
{
  init();
}


void GridMeterLocalDensity::init()
{
  m_properties.setClassName("LocalDensity");

  m_properties.setDescription("Measure the mean local density of the particles.");

//  STRINGPC
//      (weightingFunction, m_weighting_function,
//       "Weighting function to be used for local density calculation.");

  STRINGPC
      (densitySymbol, m_rhoSymbol,
       "Density to be measured.");

  BOOLPC
      (oneProp, m_oneProp,
       "Will the local density be computed over all ColourPairs (CP) or only over"
           " the CP corresponding to the chosen species?");

//  m_weighting_function = "default";
  m_oneProp = true;
  m_rhoSymbol = "undefined";
}


void GridMeterLocalDensity::setup()
{
  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterLocalDensity::setup", "Please specify a species.");

//   m_cp = M_MANAGER->cp(m_colour, m_colour);

//   m_wf = M_SIMULATION->findWeightingFunction(m_weighting_function);

  m_density_m_o = 
    M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + /*STR_LOCAL_DENSITY*/m_rhoSymbol + STR_DELIMITER + /*m_weighting_function + STR_DELIMITER +*/ "mean", DataFormat::VECTOR_DOUBLE).offset;
  m_density_v_o =
    M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + /*STR_LOCAL_DENSITY*/m_rhoSymbol + STR_DELIMITER + /*m_weighting_function + STR_DELIMITER +*/ "variance", DataFormat::VECTOR_DOUBLE).offset;

//   m_cp->registerCalc(m_density_o, new ValCalculatorRho(m_wf), m_oneProp);
  
  if(!Particle::s_tag_format[m_colour].attrExists(m_rhoSymbol))
    throw gError("GridMeterLocalDensity::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + m_species);
  m_density_o = Particle::s_tag_format[m_colour].offsetByName(m_rhoSymbol);
}

void GridMeterLocalDensity::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *values_density = data->vectorDoubleByOffset(m_density_m_o).value();
  vector<double> *variances_density = data->vectorDoubleByOffset(m_density_v_o).value();
  int n_cells = M_GRID_AVERAGER->nCells();

  values_density->resize(n_cells);
  variances_density->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       double d = __iSLFE->tag.doubleByOffset(m_density_o);

       /* fixme!!! Slow! Better divide only once per cell. */
       (*values_density)[index] += d/n;
       (*variances_density)[index] += d*d/n;
     }
    );
}

