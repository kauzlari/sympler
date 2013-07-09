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


#include "apply_vector_field.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<ApplyVectorField> apply_vector_field("ApplyVectorField");

/* --- ApplyVectorField --- */

ApplyVectorField::ApplyVectorField(Simulation* sim)
  : Thermostat(sim), m_colour(11111111)
{
  init();
}


void ApplyVectorField::init()
{
  m_properties.setClassName("ApplyVectorField");

  m_properties.setDescription(
    "When called, this callable adds a user-specified vector-field to each particle."
  );

  FUNCTIONFIXEDPC(vx, m_vecX, 
                  "This sets the x-component of the additional vector to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z'.");
  
  FUNCTIONFIXEDPC(vy, m_vecY, 
                  "This sets the y-component of the additional vector to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z'.");
  
  FUNCTIONFIXEDPC(vz, m_vecZ, 
                  "This sets the z-component of the additional vector to the specified algebraic expression. You may use constants or the known variables 'x', 'y', 'z'.");
  
  m_vecX.addVariable("x");
  m_vecX.addVariable("y");
  m_vecX.addVariable("z");
  
  m_vecY.addVariable("x");
  m_vecY.addVariable("y");
  m_vecY.addVariable("z");
  
  m_vecZ.addVariable("x");
  m_vecZ.addVariable("y");
  m_vecZ.addVariable("z");
  
  m_vecX.setExpression("0");
  m_vecY.setExpression("0");
  m_vecZ.setExpression("0");


  STRINGPC
    (species, m_species,
     "Species this callable should work on.");

  m_species = "UNDEF";

  STRINGPC
    (symbol, m_symbolName,
     "Name of the symbol to be modified.");

  m_symbolName = "undefined";

}


void ApplyVectorField::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ApplyVectorField::setup", "Attribute 'species' was not defined!");

  m_colour = M_MANAGER->getColour(m_species);

  if(m_symbolName == "undefined")
    throw gError("ApplyVectorField::setup", "Attribute 'symbol' was not defined!");

  // the attribute should already exist; the following should also search and check the format
  try
    {
      m_offset = Particle::s_tag_format[m_colour].indexOf(m_symbolName, DataFormat::POINT);
      m_offset = Particle::s_tag_format[m_colour].offsetByIndex(m_offset);
    }
  catch(gError& err)
    {
      throw gError("ApplyVectorField::setup", "search for symbol failed. The message was " + err.message()); 
    }
  
}

void ApplyVectorField::thermalize(Phase* phase)
{

    FOR_EACH_FREE_PARTICLE_C
      (phase,
       m_colour,
       point_t& vec = __iSLFE->tag.pointByOffset(m_offset);
       point_t& r = __iSLFE->r;
       vec.x += m_vecX(r.x, r.y, r.z);
       vec.y += m_vecY(r.x, r.y, r.z);
       vec.z += m_vecZ(r.x, r.y, r.z);
    );

}


