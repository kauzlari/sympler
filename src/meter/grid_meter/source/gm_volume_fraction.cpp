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



#include "gm_volume_fraction.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "grid_averager.h"
#include "val_calculator_rho.h"

GridMeter_Register<GridMeterVolumeFraction> grid_meter_volume_fraction("VolumeFraction");


#define STR_VOLUME_FRACTION string("volume_fraction")

#define M_AVERAGER  ((GridAverager*) m_parent)
#define M_SIMULATION ((Simulation*) M_AVERAGER->parent())
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

/*--- GridMeterVolumeFraction ---*/

GridMeterVolumeFraction::GridMeterVolumeFraction(GridAverager *averager): GridMeter(averager)
{
  init();
}


void GridMeterVolumeFraction::init()
{
  m_properties.setClassName("VolumeFraction");

  m_properties.setDescription(
    "Measures the density of a scalar field, i.e. phi/volume where phi "
    "is the value of the scalar field and volume is the volume per particle. The volume per particle "
    "is given by 1/d, where d is the local density. The attribute 'cutoff' defines the cut-off for the "
    "calculation of the local density."
  );

  STRINGPC
    (scalar, m_scalar_name,
     "Name of the scalar field to measure.");

  STRINGPC
      (densitySymbol, m_rhoSymbol,
       "Density to be measured.");
  
//   STRINGPC
//       (weightingFunction, m_weighting_function,
//        "Weighting function to be used for local density calculation.");

  BOOLPC
      (oneProp, m_oneProp,
       "Will the local density be computed over all ColourPairs (CP) or only over"
           " the CP corresponding to the chosen species?");

  m_scalar_name = "scalar";
//   m_weighting_function = "default";
  m_oneProp = true;
  m_rhoSymbol = "undefined";
}


void GridMeterVolumeFraction::setup()
{
  pair<size_t, size_t> temp;

  GridMeter::setup();

  if (m_colour == ALL_COLOURS)
    throw gError("GridMeterVolumeFraction::setup", "Please specify a species.");

//   m_wf = M_SIMULATION->findWeightingFunction(m_weighting_function);

  m_scalar_m_o = 
    M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + m_scalar_name + STR_DELIMITER + STR_VOLUME_FRACTION +
      /*STR_DELIMITER + m_weighting_function +*/ STR_DELIMITER + "mean", DataFormat::VECTOR_DOUBLE).offset;
  m_scalar_v_o =
    M_GRID_AVERAGER->format().addAttribute
      (m_species + STR_DELIMITER + m_scalar_name + STR_DELIMITER + STR_VOLUME_FRACTION + 
      /*STR_DELIMITER + m_weighting_function +*/ STR_DELIMITER + "variance", DataFormat::VECTOR_DOUBLE).offset;

  if(!Particle::s_tag_format[m_colour].attrExists(m_scalar_name))
    throw gError("GridMeterVolumeFraction::setup", "Symbol '" + m_scalar_name + "' not found for species " + m_species);
  m_scalar_o = Particle::s_tag_format[m_colour].offsetByName(m_scalar_name);

//   M_MANAGER->cp(m_colour, m_colour)->registerCalc(temp, new ValCalculatorRho(m_wf), m_oneProp);
//   m_density_o = temp.first; /* first and second are the same! */

  if(!Particle::s_tag_format[m_colour].attrExists(m_rhoSymbol))
    throw gError("GridMeterVolumeFraction::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + m_species);
  m_density_o = Particle::s_tag_format[m_colour].offsetByName(m_rhoSymbol);

}


void GridMeterVolumeFraction::measure(data_sp data)
{
  Phase *phase = ((Simulation*) M_GRID_AVERAGER->parent())->phase();
  vector<double> *scalar = data->vectorDoubleByOffset(m_scalar_m_o).value();
  vector<double> *variances = data->vectorDoubleByOffset(m_scalar_v_o).value();
  int n_cells = M_GRID_AVERAGER->nCells();

  scalar->resize(n_cells);
  variances->resize(n_cells);

  FOR_EACH_FREE_PARTICLE_C
    (phase, m_colour,
     int index = M_GRID_AVERAGER->location(c, __iSLFE->mySlot);
     int n = M_GRID_AVERAGER->nParticles(c, index);
         
     if (n) {
       double s = __iSLFE->tag.doubleByOffset(m_scalar_o)*__iSLFE->tag.doubleByOffset(m_density_o);

       /* fixme!!! Slow! Better divide only once per cell. */
       (*scalar)[index] += s/n;
       (*variances)[index] += s*s/n;
     }
    );
}

