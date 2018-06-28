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



#include "apply_vel_field.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<ApplyVelField> apply_vel_field("ApplyVelField");

/* --- ApplyVelField --- */

ApplyVelField::ApplyVelField(Simulation* sim)
  : Thermostat(sim), m_colour(11111111)
{
  init();
}


void ApplyVelField::init()
{
  m_properties.setClassName("ApplyVelField");

  m_properties.setDescription(
    "When called, this callable modifies the velocity-field of each particle according to user-defined expressions for the attributes 'u', 'v', 'w' described below. In the expressions you may use constants, functions (listed under sympler --help expressions) and the known variables 'x', 'y', 'z', 'u', 'v', 'w' (NOT those listed under sympler --help expressions!). \nNOTE: The velocity is SET and not incremented. If you want to increment, use the known variables 'u', 'v', 'w' containing the respective old values of the velocity components, e.g., by typing 'u = \"u + exp(x)\"', where the second part is your increment." 
  );

  FUNCTIONFIXEDPC(u, m_velX, 
                  "This sets the x-component of the velocity to the specified algebraic expression. See above for allowed expressions.");
  
  FUNCTIONFIXEDPC(v, m_velY, 
                  "This sets the y-component of the velocity to the specified algebraic expression. See above for allowed expressions.");
  
  FUNCTIONFIXEDPC(w, m_velZ, 
                  "This sets the z-component of the velocity to the specified algebraic expression. See above for allowed expressions.");
  
  m_velX.addVariable("x");
  m_velX.addVariable("y");
  m_velX.addVariable("z");
  
  m_velY.addVariable("x");
  m_velY.addVariable("y");
  m_velY.addVariable("z");
  
  m_velZ.addVariable("x");
  m_velZ.addVariable("y");
  m_velZ.addVariable("z");

  m_velX.addVariable("u");
  m_velX.addVariable("v");
  m_velX.addVariable("w");
  
  m_velY.addVariable("u");
  m_velY.addVariable("v");
  m_velY.addVariable("w");
  
  m_velZ.addVariable("u");
  m_velZ.addVariable("v");
  m_velZ.addVariable("w");
  
  m_velX.setExpression("u");
  m_velY.setExpression("v");
  m_velZ.setExpression("w");


  STRINGPC
    (species, m_species,
     "Species this callable should work on.");

  m_species = "UNDEF";

}


void ApplyVelField::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ApplyVelField::setup", "Attribute 'species' was not defined!");

  m_colour = M_MANAGER->getColour(m_species);
}

void ApplyVelField::thermalize(Phase* phase)
{

    FOR_EACH_FREE_PARTICLE_C
      (phase,
       m_colour,
       __iSLFE->v.x = m_velX(__iSLFE->r.x, __iSLFE->r.y, __iSLFE->r.z, __iSLFE->v.x, __iSLFE->v.y, __iSLFE->v.z);
       __iSLFE->v.y = m_velY(__iSLFE->r.x, __iSLFE->r.y, __iSLFE->r.z, __iSLFE->v.x, __iSLFE->v.y, __iSLFE->v.z);
       __iSLFE->v.z = m_velZ(__iSLFE->r.x, __iSLFE->r.y, __iSLFE->r.z, __iSLFE->v.x, __iSLFE->v.y, __iSLFE->v.z);
    );

}


