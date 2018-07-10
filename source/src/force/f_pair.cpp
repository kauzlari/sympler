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



#include "manager_cell.h"
#include "simulation.h"
#include "weighting_function.h"
#include "colour_pair.h"

#include "f_pair.h"


#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

FPair::FPair(Simulation *simulation): GenF(simulation)/*, m_cp(NULL)*/
{
  init();
}


FPair::~FPair()
{
}


void FPair::init()
{
  m_properties.setName("FPair");

  m_properties.setDescription("Base class for all pair forces.");
    
  STRINGPC
    (species1, m_species.first,
     "First species, this force should act on.");

  STRINGPC
    (species2, m_species.second,
     "Second species, this force should act on.");

  m_species.first = "UNDEF";
  m_species.second = "UNDEF";

  m_is_pair_force = true;
  m_is_particle_force = false;

}


void FPair::setup()
{
  GenF::setup();
  
  if (m_species.first == "UNDEF")
    throw gError("FPairVector::setup", "'species1' undefined in " + m_properties.name());
  if (m_species.second == "UNDEF")
    throw gError("FPairVector::setup", "'species2' undefined in " + m_properties.name());

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  
  m_cp->registerForce(this);

  m_cp->setNeedPairs(true);
	
}

