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


#include "val_calculator_dirichlet_BC_vels.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"

// #include <utility>

const SymbolRegister<ValCalculatorDirichletBCVels> val_calc_dirichlet_BC_vels("DirichletBCVels");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorDirichletBCVels::ValCalculatorDirichletBCVels(string symbol)
  : ValCalculatorBC(symbol)
{
  m_stage = 0;
  //   MSG_DEBUG("ValCalculatorDirichletBCVels::ValCalculatorDirichletBCVels", "CONSTRUCTOR");
}

ValCalculatorDirichletBCVels::ValCalculatorDirichletBCVels(/*Node*/Simulation* parent)
  : ValCalculatorBC(parent)
{
  m_stage = 0;
  init();
}

void ValCalculatorDirichletBCVels::init()
{
  m_properties.setClassName("DirichletBCVels");

  m_properties.setDescription("Saves the pair-specific value of the velocity of the boundary particle used for applying a Dirichlet boundary condition (BC) in each pair of particles. This calculator uses a linear approximation. The actual value of the Dirichlet boundary condition is assumed to be stored in the velocity field of the respective boundary particle.");
  
//   STRINGPC
//       (symbol, m_symbolName,
//        "Name of the symbol for the boundary value.");
  
//   STRINGPC(wallSpecies, m_wallSpecies, 
//            "Species of the wall particles."
//           );

  m_symbolName = "vBC";
  m_wallSpecies = "undefined";
  
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorDirichletBCVels::setup()
{
//   M_CONTROLLER->registerForSetupAfterParticleCreation(this);
//   if(m_species.first == "undefined")
//     throw gError("ValCalculatorDirichletBCVels::setup", "Attribute 'species1' has value \"undefined\" .");
//   if(m_species.second == "undefined")
//     throw gError("ValCalculatorDirichletBCVels::setup", "Attribute 'species2' has value \"undefined\" .");
//   ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

//   // check whether wall particles are first or second in the pairs
//   if(m_wallSpecies == cp->firstSpecies())
//     m_wallIsSecond = false;
//   else if(m_wallSpecies == cp->secondSpecies())
//     m_wallIsSecond = true;
//   else
//     throw gError("ValCalculatorDirichletBCVels::setup", "Attribute wallSpecies = \"" + m_wallSpecies + "\" is neither equal to attribute species1 = \"" + m_species.first + "\" nor to species2 = \"" + m_species.second + "\".");

//   m_wallColour = M_MANAGER->getColour(m_wallSpecies);

  m_datatype = DataFormat::POINT;
  
  ValCalculatorBC::setup();
  
}

// void /*pair<size_t, size_t>*/ ValCalculatorDirichletBCVels::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
// {
//   MSG_DEBUG("ValCalculatorDirichletBCVels::setSlot", "CALLED");
//   m_slot = slot = cp->tagFormat().addAttribute
//       ("ValCalculator_" + myName() + "_" + cp->toString(), DataFormat::POINT, false, m_symbolName).offset;
// }


