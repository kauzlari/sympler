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



#include "lennard_jones_vc.h"

#include "phase.h"
#include "simulation.h"
#include "val_calculator_r_i6.h"
#include "val_calculator_r_i.h"
#include "pair_creator.h"
#include "threads.h"


using namespace std;

const GenFTypeConcr<LJVC> lj_vc("LJVC");


#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()


//---- Constructors/Destructor ----

LJVC::LJVC(Simulation *simulation)
    : LJ(simulation)
{
    init();
}

void LJVC::init()
{
  m_properties.setClassName("LJVC");

  m_properties.setDescription(
    "Force based on the Lennard-Jones potential. Parametrization is "
    "chosen to be: V(r)=4*epsilon*[(sigma/r)^12-(sigma/r)^6]. This modules "
    "creates two Symbol modules storing 1/r and 1/r^6 for each pair, which "
    "may be used by other modules."
  );
}


void LJVC::setup()
{
  LJ::setup();
  ColourPair *cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  cp->registerCalc(m_compute_r6i_offset, new ValCalculatorRi6, false);
  cp->registerCalc(m_compute_ri_offset, new ValCalculatorRi, false);

  m_c1 = M_MANAGER->getColour(m_species.first);
  m_c2 = M_MANAGER->getColour(m_species.second);

}

void LJVC::computeForces(Particle* part, int force_index)
{
  throw gError("LJVC::computeForces", "Fatal error: do not call FPairVels::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void LJVC::computeForces(int force_index)
{
  throw gError("LJVC::computeForces", "Fatal error: do not call FPairVels::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}
