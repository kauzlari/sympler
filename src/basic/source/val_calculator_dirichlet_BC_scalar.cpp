/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


#include "val_calculator_dirichlet_BC_scalar.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"
#include "particle_cache.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"

// #include <utility>

const SymbolRegister<ValCalculatorDirichletBCScalar> val_calc_dirichlet_BC_scalar("DirichletBCScalar");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorDirichletBCScalar::ValCalculatorDirichletBCScalar(string symbol)
  : ValCalculatorBC(symbol)
{
//   MSG_DEBUG("ValCalculatorDirichletBCScalar::ValCalculatorDirichletBCScalar", "CONSTRUCTOR");
}

ValCalculatorDirichletBCScalar::ValCalculatorDirichletBCScalar(/*Node*/Simulation* parent)
  : ValCalculatorBC(parent)
{
  init();
}

void ValCalculatorDirichletBCScalar::init()
{
  m_properties.setClassName("DirichletBCScalar");

  m_properties.setDescription("Saves the pair-specific value of an arbitrary scalar of the boundary particle used for applying a Dirichlet boundary condition (BC) in each pair of particles. This calculator uses a linear approximation as described in [70]. The actual value of the Dirichlet boundary condition is assumed to be stored in the symbol of the respective boundary particle given by the attribute 'scalar'.\n"
"[70]: J. P. Morris, P. J. Fox, Y. Zhu, J. Comp. Phys. 136 (1997) 214–226.");
  
   STRINGPC
       (scalar, m_scalarName,
        "Name of the scalar the BC should be applied to.");
  
  
//   STRINGPC(wallSpecies, m_wallSpecies, 
//            "Species of the wall particles."
//           );

   m_symbolName = "undefined";
   m_scalarName = "undefined";
   m_wallSpecies = "undefined";
  
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorDirichletBCScalar::setup()
{
  if(m_symbolName == "undefined")
    throw gError("ValCalculatorDirichletBCScalar::setup", "Attribute 'symbol' has value \"undefined\" .");

  if(m_scalarName == "undefined")
    throw gError("ValCalculatorDirichletBCScalar::setup", "Attribute 'scalar' has value \"undefined\" .");

  m_datatype = DataFormat::DOUBLE;
  
  ValCalculatorBC::setup();  

  // now (2016-12-14), ValCalculatorBC::setup() has defined m_species
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  try {
    DataFormat::attribute_t firstAttr =
      Particle::s_tag_format[cp->firstColour()].attrByName(m_scalarName);
  if(firstAttr.datatype != DataFormat::DOUBLE)
    throw gError(/*"ValCalculatorDirichletBCScalar::setup"*/"", "the symbol " + m_scalarName +
		 " is registerd as a non-scalar for species " +
		 cp->manager()->species(cp->firstColour()));
  }
  catch(gError& err) {
    throw gError("ValCalculatorDirichletBCScalar::setup", "For species " + cp->firstSpecies() + ": " + err.msg());
  }

  try {
    DataFormat::attribute_t secondAttr =
      Particle::s_tag_format[cp->secondColour()].attrByName(m_scalarName);
  if(secondAttr.datatype != DataFormat::DOUBLE)
    throw gError(/*"ValCalculatorDirichletBCScalar::setup"*/ "", "the symbol " + m_scalarName +
		 " is registerd as a non-scalar for species " +
		 cp->manager()->species(cp->secondColour()));
  }
  catch(gError& err) {
    throw gError("ValCalculatorDirichletBCScalar::setup", "For species " + cp->secondSpecies() + ": " + err.msg());
  }
     
   m_scalarOffset.first = 
     Particle::s_tag_format[cp->firstColour()].offsetByName(m_scalarName);

   m_scalarOffset.second = 
     Particle::s_tag_format[cp->secondColour()].offsetByName(m_scalarName);
  
}

