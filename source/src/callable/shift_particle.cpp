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

#include "shift_particle.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<ShiftParticle> shift_particle("ShiftParticle");

/* --- ShiftParticle --- */

ShiftParticle::ShiftParticle(Simulation* sim)
  : Thermostat(sim), m_colour(11111111) 
{
  init();
}


void ShiftParticle::init()
{
  m_properties.setClassName("ShiftParticle");

  m_properties.setDescription(
    "Callable which shifts particles by a user-defined particle-"
    "specific displacement vector specified by the symbol given in "
    "attribute 'shiftSymbol'."
    // "\nNOTE1: Do not define any Callables after (below) ShiftParticle, "
    // "which require an up-to-date neighbour list, because this "
    // "ShiftParticle changes particle positions but DOES NOT invoke a "
    // "neighbour list update and hence no recomputation of particle "
    // "distances."
  );

  STRINGPC
    (species, m_species,
     "Species this callable should work on.");

  m_species = "UNDEF";

  STRINGPC
    (shiftSymbol, m_shiftSymbolName,
     "Name of the vector symbol holding the shift for each particle.");

  m_shiftSymbolName = "undefined";

}


void ShiftParticle::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ShiftParticle::setup", "Attribute 'species' was not "
		 "defined!");

  m_colour = M_MANAGER->getColour(m_species);

  if(m_shiftSymbolName == "undefined")
    throw gError("ShiftParticle::setup", "Attribute 'shiftSymbol' has "
		 "value \"undefined\"");

    // the attribute should already exist
    try
    {
      m_shiftOffset =
	Particle::s_tag_format[m_colour].indexOf
	(m_shiftSymbolName, DataFormat::POINT);
      
      m_shiftOffset =
	Particle::s_tag_format[m_colour].offsetByIndex(m_shiftOffset);
    }
    catch(gError& err)
    {
      throw gError("ShiftParticle::setup", "search for symbol given "
		   "by attribute 'shiftSymbol' failed. The message "
		   "was: " + err.message()); 
    }

}

void ShiftParticle::thermalize(Phase* phase)
{

  Controller* controller = M_CONTROLLER;
  
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase,
     m_colour,
     this,
     
     i->r += i->tag.pointByOffset(m_shiftOffset);
     
     );

  // notify cells and neighbour list of new positions
  // FIXME: Currently we do not pass the commented out argument, which
  // would be a way to let the Phase and the Cells check for wall
  // collisions based on the shift. If there is a hit, the walls could
  // move the particle back along the line r_vec - shift_vec until it
  // is back inside and a "sufficiently large epsilon" away from the
  // walls. This is not yet (2018-04-17) implemented within the
  // Cell-Wall logic
  phase -> invalidatePositions(m_colour/*, m_shiftOffset*/);

  // update neighbour-list such that shifting immediately effective
  // for further computations requiring neighbour information
  controller -> triggerNeighbourUpdate();
  
}


