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
  #include <freesteam/b23.h>
  #include <freesteam/steam_pT.h>
  #include <freesteam/region3.h>
}

#include "pca_iapws-if97_p.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<PCacheIAPWSIF97p> PCacheIAPWSIF97p("PCacheIAPWSIF97p");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCacheIAPWSIF97p::PCacheIAPWSIF97p
  (size_t colour, size_t offset, string symbolName)
  : PCacheIAPWSIF97TwoVar(colour, offset, symbolName) {
}

PCacheIAPWSIF97p::PCacheIAPWSIF97p
    (/*Node*/Simulation* parent)
  : PCacheIAPWSIF97TwoVar(parent) {
  init();
}


void PCacheIAPWSIF97p::checkConstraints() {

  if (m_var2Min <= 623.15 || m_var2Max >= 863.15)
     throw gError
       ("PCacheIAPWSIF97p::setup", "Requested temperature parameters "
	"aren't within region 3 (623.15K < T < 863.15K). (See "
	"IAPWS-IF97 for more information)");

}


void PCacheIAPWSIF97p::freesteamCalculationForState(double& result) const {

  // Calculates the minimum and maximum density boarders to check
  // if the given density is in Region 3.
  double& Tin = *(m_inputVarPtrs[1]);
  double& densityIn = *(m_inputVarPtrs[0]);
  
  double b23Pressure = freesteam_b23_p_T(Tin);
  SteamState S = freesteam_set_pT(b23Pressure, Tin);
  double densityBoundary = freesteam_rho(S);
  
  if (densityIn > densityBoundary) {
    result =
      freesteam_region3_p_rhoT(densityIn, Tin/*temperature*/);
  }
  else {
    throw gError("PCacheIAPWSIF97::setup", "Requested density and temperature parameters aren't within region 3. (See IAPWS-IF97 for more information): Unsuccessfully tried to compute P(rho=" + ObjToString(*(m_inputVarPtrs[0])) + ",T=" + ObjToString(Tin) + ") crossing rho_min(T) = " + ObjToString(densityBoundary));
  }
  
}


void PCacheIAPWSIF97p::init()
{
  m_properties.setClassName("PCacheIAPWSIF97p");
  m_properties.setName("PCacheIAPWSIF97p");
  m_properties.setDescription
    (m_properties.description() +
     "\nFor module PCacheIAPWSIF97p, var1 = density (rho), var2 = "
     "temperature (T), and out = pressure (p)."
     "NOTE: This module can currently only compute pressure, if "
     "(rho,T) is within region 3. See IAPWS-IF97 for more "
     "information."
     ); 
}


void PCacheIAPWSIF97p::setup()
{
  PCacheIAPWSIF97TwoVar::setup();
}


void PCacheIAPWSIF97p::registerWithParticle()
{
}

#endif    // HAVE_FREESTEAM

