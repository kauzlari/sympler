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



#include "integrator.h"

#include "phase.h"
#include "simulation.h"
#include "controller.h"
#include "manager_cell.h"


REGISTER_SMART_ENUM
(Integrator_Factory,
 "Integrators advance degrees of freedom for every timestep."
);

#define M_CONTROLLER  ((Controller*) m_parent)
#define M_SIMULATION  ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

Integrator::Integrator(Controller *controller): Node((Node*) controller)
{
  init();
}


Integrator::~Integrator()
{
}


void Integrator::init()
{
	STRINGPC
    (species, m_species,
     "Species, this integrator is intended for.");
	
  m_species = "UNDEF";
#ifdef _OPENMP
  m_vec_offset.resize(global::n_threads);
  m_merge = false;
#endif
}

size_t Integrator::getColourAndAdd(string species) {
	return M_MANAGER->getColourAndAdd(species);
}

void Integrator::setup()
{
  Node::setup();

  if (m_species == "UNDEF") {
    throw gError("Integrator::setup", "No setup for species UNDEF, please define species.");
  }

  m_colour = M_MANAGER->getColourAndAdd(m_species);

  // next is needed because setClassName sets both name and className
  string name = m_properties.name();

  m_properties.setClassName(m_properties.className() + "(" + m_species + ")");
  m_properties.setName(name + "(" + m_species + ")");

  MSG_DEBUG("Integrator::setup", "m_colour = " << m_colour);
}
