/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


extern "C" {
  #include <freesteam/steam_pT.h>
}
#include "density_calculation.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<DensityCalculation> DensityCalculation("DensityCalculation");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


DensityCalculation::DensityCalculation
  (size_t colour, size_t offset, string symbolName)
  : PCacheIAPWSIF97(colour, offset, symbolName) {
}

DensityCalculation::DensityCalculation
    (/*Node*/Simulation* parent)
  : PCacheIAPWSIF97(parent) {
  init();
}


void DensityCalculation::checkConstraints() {}


void DensityCalculation::freesteamCalculationForState
(double& result, const double& inputVar1, const double& inputVar2)
  const {
  
  SteamState S = freesteam_set_pT(inputVar1, inputVar2);
  result = freesteam_rho(S);
  
}


void DensityCalculation::init()
{
  m_properties.setClassName("DensityCalculation");
  m_properties.setName("DensityCalculation");
  m_properties.setDescription
    (m_properties.description() +
     "\nFor module DensityCalculation, var1 = pressure (p), var2 = "
     "temperature (T), and out = density (rho)."
     ); 
}


void DensityCalculation::setup()
{
  PCacheIAPWSIF97::setup();
}


void DensityCalculation::registerWithParticle() {}
