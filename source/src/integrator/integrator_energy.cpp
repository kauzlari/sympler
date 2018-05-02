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



#include "gen_f.h"
#include "phase.h"
#include "particle.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_energy.h"
#include "threads.h"

using namespace std;


#define M_CONTROLLER  ((Controller*) m_parent)
#define M_SIMULATION  ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()


const Integrator_Register<IntegratorEnergy> integrator_energy("IntegratorEnergy");

//---- Constructors/Destructor ----

IntegratorEnergy::IntegratorEnergy(Controller *controller): IntegratorScalar(controller)
{
  init();
}


IntegratorEnergy::~IntegratorEnergy()
{
}



//---- Methods ----

void IntegratorEnergy::init()
{
//   MSG_DEBUG("IntegratorEnergy::init()", "running");
  m_properties.setClassName("IntegratorEnergy");

  m_properties.setDescription(
    "Adds an additional degree of freedom, the internal energy, "
    "to the particles specified. This allows for non-isothermal DPD."
      " Integration is performed with a simple Euler scheme."
  );

  /* Register the energy degree of freedom for the particles. */

  FUNCTIONFIXEDPC(ds_de, m_ds_de,
    "Differential of the entropy with respect to the energy."
    " Here, the variable 'e' is allowed and represents the internal energy.\n"
    "Example: ds_de = '\"1000/e\"\nThis would mean that "
    "e = 1000/(ds_de) = 1000*T, where T = (ds_de)^-1 is the temperature. "
    "This implies a constant heat capacity of 1000.");

  m_ds_de.addVariable("e");

  m_scalar_name = "internal_energy";
  m_scalar_symbol = "e";
}

void IntegratorEnergy::setup()
{
  IntegratorScalar::setup();

  m_T_offset =
    Particle::s_tag_format[m_colour].addAttribute
      ("temperature",
       DataFormat::DOUBLE,
       true,
       "T").offset;
  MSG_DEBUG("IntegratorEnergy::setup", "m_colour = " << m_colour << ", T_offset = " << m_T_offset);
}

// calculate the temperature
void IntegratorEnergy::deriveQuantities()
{
/*  MSG_DEBUG("IntegratorEnergy::deriveQuantities", "m_colour = " << m_colour << ", T_offset = " << m_T_offset);*/
  Phase *phase = M_PHASE;

  // computation of temperature
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     i->tag.doubleByOffset(((IntegratorEnergy*) data)->m_T_offset) =
     1/((IntegratorEnergy*) data)->dEntropy_dEnergy
     (
       i->tag.doubleByOffset(((IntegratorEnergy*) data)->m_scalar_offset)
     );
     assert(i->tag.doubleByOffset(((IntegratorEnergy*) data)->m_T_offset) != 0);
    );
}


#ifdef _OPENMP
string IntegratorEnergy::dofIntegr() {
  return "temperature";
}


#endif
