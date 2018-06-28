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


#include "apply_scalar_field.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<ApplyScalarField> apply_scalar_field("ApplyScalarField");

/* --- ApplyScalarField --- */

ApplyScalarField::ApplyScalarField(Simulation* sim)
  : Thermostat(sim), m_colour(11111111)
{
  init();
}


void ApplyScalarField::init()
{
  m_properties.setClassName("ApplyScalarField");

  m_properties.setDescription(
    "When called, this callable modifies the user-specified scalar-field of each particle according to a user-defined expression for the attributes 'scalar' described below. In the expressions you may use constants, functions (listed under sympler --help expressions) and the known variables 'x', 'y', 'z', 'scalar' (NOT those listed under sympler --help expressions!). \nNOTE: The scalar is SET and not incremented. If you want to increment, use the known variable 'scalar' containing the old value, e.g., by typing 'scalar = \"scalar + exp(x)\"', where the second part is your increment." 
  );

  FUNCTIONFIXEDPC(scalar, m_expr, 
                  "This sets the scalar field to the algebraic expression specified here. See above for allowed expressions. Note that the affected scalar-field itself is accessed by the variable-name 'scalar' and NOT by its original name given in attribute 'symbol'.");
  
  m_expr.addVariable("x");
  m_expr.addVariable("y");
  m_expr.addVariable("z");
  
  m_expr.addVariable("scalar");

  // scalar unaltered by default
  m_expr.setExpression("scalar");
 
  STRINGPC
    (species, m_species,
     "Species this callable should work on.");

  m_species = "UNDEF";

  STRINGPC
    (symbol, m_symbolName,
     "Name of the symbol to be modified.");

  m_symbolName = "undefined";

}

void ApplyScalarField::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ApplyScalarField::setup", "Attribute 'species' was not defined!");

  m_colour = M_MANAGER->getColour(m_species);

  if(m_symbolName == "undefined")
    throw gError("ApplyScalarField::setup", "Attribute 'symbol' was not defined!");

  // the attribute should already exist; the following should also search and check the format
  try {
    m_offset = Particle::s_tag_format[m_colour].indexOf(m_symbolName, DataFormat::DOUBLE);
    m_offset = Particle::s_tag_format[m_colour].offsetByIndex(m_offset);
  }
  catch(gError& err) {
    throw gError("ApplyScalarField::setup", "search for symbol failed. The message was " + err.message()); 
  }
  
}

void ApplyScalarField::thermalize(Phase* phase)
{
  FOR_EACH_FREE_PARTICLE_C
    (phase,
     m_colour,
     double& scalar = __iSLFE->tag.doubleByOffset(m_offset);
     point_t& r = __iSLFE->r;
     
     scalar = m_expr(r.x, r.y, r.z, scalar);
     );
}


