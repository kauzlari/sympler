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



#include <algorithm>

#include "pair_creator.h"
#include "threads.h"
#include "inlet_cell.h"
#include "simulation.h"
#include "vertex_list.h"
#include "manager_cell.h"
#include "particle_cache.h"
#include "particle_creator.h"
#include "colour_pair.h"

using namespace std;

int PairCreator::counterTN = 0;



REGISTER_SMART_ENUM
(PairCreator_Factory,
 "As the name suggests, PairCreators define the way the particle pairs should be created. "
 "They are children of a phase."
);


//---- Constructors/Destructor ----

PairCreator::PairCreator(Phase* p)
  : Node((Node*) p), m_valid_dist(false)
{
  init();
}


PairCreator::~PairCreator()
{
}


void PairCreator::init()
{
  m_properties.setClassName("PairCreator");

  m_properties.setDescription("Determines the way the pairDistances are created.");
}


void PairCreator::setup()
{
  Node::setup();
}

void PairCreator::createDistances()
{
}

void PairCreator::invalidatePositions()
{
}

void PairCreator::pairAddedFree(Pairdist* pd)
{
}

void PairCreator::pairAddedFrozen(Pairdist* pd)
{
}


