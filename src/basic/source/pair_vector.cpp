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



#include "pair_vector.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

const SymbolRegister<PairVector> pair_vector("PairVector");

PairVector::PairVector(string symbol) :
  PairArbitrary(symbol) {
}
////Constructor/Destructor
PairVector::PairVector(Simulation* parent) :
  PairArbitrary(parent) {
  m_function.setReturnType(Variant::VECTOR);
  
  
  m_datatype = DataFormat::POINT;
  init();
}

PairVector::~PairVector() {
}

void PairVector::init() {
  m_properties.setClassName("PairVector");
  m_properties.setDescription("Module for calculation of the pair factor between a pair of particles "
			      "saves it in the pair attribute");
  
}

//---- Methods ----

void PairVector::setup() {
  
  MSG_DEBUG("PairVector::setup()","Set up of PairVector");
  PairArbitrary::setup();
  
}

