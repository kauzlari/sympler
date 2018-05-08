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


#ifdef HAVE_FREESTEAM

extern "C" {
  #include <freesteam/viscosity.h>
}

#include "pca_iapws-if97_eta.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<PCacheIAPWSIF97eta> PCacheIAPWSIF97eta("PCacheIAPWSIF97eta");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCacheIAPWSIF97eta::PCacheIAPWSIF97eta
  (size_t colour, size_t offset, string symbolName)
  : PCacheIAPWSIF97TwoVar(colour, offset, symbolName) {
}

PCacheIAPWSIF97eta::PCacheIAPWSIF97eta
    (/*Node*/Simulation* parent)
  : PCacheIAPWSIF97TwoVar(parent) {

   init();
}


void PCacheIAPWSIF97eta::checkConstraints() {

}


void PCacheIAPWSIF97eta::freesteamCalculationForState(double& result) const {

  // assuming that m_inputVarPtrs[0..1] were set correctly before
  // calling this function
  result =
    freesteam_mu_rhoT(*(m_inputVarPtrs[0])/*density*/,
		      *(m_inputVarPtrs[1])/*temperature*/);

}


void PCacheIAPWSIF97eta::init()
{
  m_properties.setClassName("PCacheIAPWSIF97eta");
  m_properties.setName("PCacheIAPWSIF97eta");
  m_properties.setDescription
    (m_properties.description() +
     "\nFor module PCacheIAPWSIF97eta, var1 = density (rho), var2 = "
     "temperature (T), and out = dynamic viscosity (eta)."
     ); 
}


void PCacheIAPWSIF97eta::setup()
{
  PCacheIAPWSIF97TwoVar::setup();
}


void PCacheIAPWSIF97eta::registerWithParticle()
{
}

#endif    // HAVE_FREESTEAM

