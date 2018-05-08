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
  #include <freesteam/steam_pT.h>
  #include <freesteam/derivs.h>
}

#include "pca_iapws-if97_cp.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<PCacheIAPWSIF97Cp> PCacheIAPWSIF97Cp("PCacheIAPWSIF97Cp");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


// char[3] s_xyz = {'T', 'p', 'h'} is used for computing
// dz/dx|y = dh/dT|p
const char PCacheIAPWSIF97Cp::s_xyz[3] = {'T', 'p', 'h'}; 


PCacheIAPWSIF97Cp::PCacheIAPWSIF97Cp
  (size_t colour, size_t offset, string symbolName)
  : PCacheIAPWSIF97OneVar(colour, offset, symbolName) {
}

PCacheIAPWSIF97Cp::PCacheIAPWSIF97Cp
    (/*Node*/Simulation* parent)
  : PCacheIAPWSIF97OneVar(parent) {
  init();
}


void PCacheIAPWSIF97Cp::checkConstraints() {}


void PCacheIAPWSIF97Cp::freesteamCalculationForState
(double& result) const {

  // assuming that m_inputVarPtrs[0] was set correctly before calling
  // this function
  SteamState S
    = freesteam_set_pT(m_constP, *(m_inputVarPtrs[0])/*temperature*/);

  // char[3] s_xyz = {'T', 'p', 'h'} is used for computing
  // dz/dx|y = dh/dT|p
  // (char*) conversion because s_xyz is const char* which I currently
  // (2018-03-29) think to be a good idea.  
  result = freesteam_deriv(S, (char*) s_xyz);
  
}


void PCacheIAPWSIF97Cp::init()
{
  m_properties.setClassName("PCacheIAPWSIF97Cp");
  m_properties.setName("PCacheIAPWSIF97Cp");
  m_properties.setDescription
    (m_properties.description() +
     "\nFor module PCacheIAPWSIF97Cp, var1 = temperature (T), "
     "and the constant parameter represents pressure (P). "
     "The output is the heat capacity at constant pressure (Cp)."
     ); 

  DOUBLEPC(constP, m_constP, 0.,
	   "Value of the constant parameter");
  
  m_constP = -HUGE_VAL;
}


void PCacheIAPWSIF97Cp::setup()
{

  PCacheIAPWSIF97OneVar::setup();

  if(m_constP == -HUGE_VAL)
    throw gError("PCacheIAPWSIF97Cp::setup", "It seems you have not defined a reasonable value for attribute 'constP'.");
  
}


void PCacheIAPWSIF97Cp::registerWithParticle() {}

#endif    // HAVE_FREESTEAM
