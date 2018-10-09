/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018,
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#include "val_calculator_dirichlet_BC_arbitrary.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"
#include "particle_cache.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorDirichletBCArbitrary::ValCalculatorDirichletBCArbitrary(string symbol)
  : ValCalculatorBC(symbol)
{
//   MSG_DEBUG("ValCalculatorDirichletBCArbitrary::ValCalculatorDirichletBCArbitrary", "CONSTRUCTOR");
}

ValCalculatorDirichletBCArbitrary::ValCalculatorDirichletBCArbitrary(/*Node*/Simulation* parent)
  : ValCalculatorBC(parent)
{
  init();
}

void ValCalculatorDirichletBCArbitrary::init()
{
	m_properties.setClassName("DirichletBCArbitrary");
  
  STRINGPC
       (variable, m_varName,
        "Name of the variable the BC should be applied to.");
  
  
//   STRINGPC(wallSpecies, m_wallSpecies, 
//            "Species of the wall particles."
//           );

   m_symbolName = "undefined";
   m_varName = "undefined";
   m_wallSpecies = "undefined";
  
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorDirichletBCArbitrary::setup()
{
  if(m_symbolName == "undefined")
    throw
			gError("ValCalculatorDirichletBCArbitrary::setup", " for module "
    		+ this->className()
				+ ": Attribute 'symbol' has value \"undefined\" .");

  if(m_varName == "undefined")
    throw
			gError("ValCalculatorDirichletBCArbitrary::setup", " for module "
					+ this->className()
					+ ": Attribute 'variable' has value \"undefined\" .");
  
  ValCalculatorBC::setup();  

  // now (2016-12-14), ValCalculatorBC::setup() has defined m_species
  ColourPair* cp
		= M_MANAGER->cp(M_MANAGER->getColour(m_species.first)
  		, M_MANAGER->getColour(m_species.second)/*m_species*/);

  try {
    DataFormat::attribute_t firstAttr =
      Particle::s_tag_format[cp->firstColour()].attrByName(m_varName);
  if(firstAttr.datatype != m_datatype)
    throw gError("ValCalculatorDirichletBCArbitrary::setup", " for module "
    		+ this->className() + ": the symbol " + m_varName
				+ " is already registered as a non "
    		+ DataFormat::attribute_t::datatypeAsString(m_datatype)
  			+ " with wrong datatype "
				+ DataFormat::attribute_t::datatypeAsString(firstAttr.datatype)
  			+ " for species " + cp->manager()->species(cp->firstColour()));
  }
  catch(gError& err) {
    throw
			gError("ValCalculatorDirichletBCArbitrary::setup", "For species "
					+ cp->firstSpecies() + ": " + err.msg());
  }

  try {
    DataFormat::attribute_t secondAttr =
      Particle::s_tag_format[cp->secondColour()].attrByName(m_varName);
  if(secondAttr.datatype != m_datatype)
    throw gError("ValCalculatorDirichletBCArbitrary::setup", " for module "
    		+ this->className() + ": the symbol " + m_varName
		    + " is already registered as a non "
    		+ DataFormat::attribute_t::datatypeAsString(m_datatype)
  			+ " with wrong datatype "
				+ DataFormat::attribute_t::datatypeAsString(secondAttr.datatype)
  			+ " for species " + cp->manager()->species(cp->secondColour()));
  }
  catch(gError& err) {
    throw
			gError("ValCalculatorDirichletBCArbitrary::setup", " for module "
					+ this->className() + ": For species " + cp->secondSpecies()
					+ ": " + err.msg());
  }
     
   m_varOffset.first =
     Particle::s_tag_format[cp->firstColour()].offsetByName(m_varName);

   m_varOffset.second =
     Particle::s_tag_format[cp->secondColour()].offsetByName(m_varName);
  
}

