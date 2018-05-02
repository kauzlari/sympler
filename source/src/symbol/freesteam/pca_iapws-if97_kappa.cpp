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
  #include <freesteam/thcond.h>
}

#include "pca_iapws-if97_kappa.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<PCacheIAPWSIF97kappa> PCacheIAPWSIF97kappa("PCacheIAPWSIF97kappa");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCacheIAPWSIF97kappa::PCacheIAPWSIF97kappa
  (size_t colour, size_t offset, string symbolName)
  : PCacheIAPWSIF97TwoVar(colour, offset, symbolName) {
}

PCacheIAPWSIF97kappa::PCacheIAPWSIF97kappa
    (/*Node*/Simulation* parent)
  : PCacheIAPWSIF97TwoVar(parent) {
  init();
}


void PCacheIAPWSIF97kappa::checkConstraints() {

}


void PCacheIAPWSIF97kappa::freesteamCalculationForState(double& result) const {

  result = freesteam_k_rhoT(*(m_inputVarPtrs[0])/*density*/,
			    *(m_inputVarPtrs[1])/*temperature*/);
}


void PCacheIAPWSIF97kappa::init()
{
  m_properties.setClassName("PCacheIAPWSIF97kappa");
  m_properties.setName("PCacheIAPWSIF97kappa");
  m_properties.setDescription
    (m_properties.description() +
     "\nFor module PCacheIAPWSIF97kappa, var1 = density (rho), var2 = "
     "temperature (T), and out = thermal conductivity (kappa)."
     ); 
}


void PCacheIAPWSIF97kappa::setup()
{
  PCacheIAPWSIF97TwoVar::setup();
}


void PCacheIAPWSIF97kappa::registerWithParticle()
{
}

