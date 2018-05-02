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


#include "bonded_pair_scalar.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

const SymbolRegister<BondedPairScalar> bonded_pair_scalar("BondedPairScalar");

BondedPairScalar::BondedPairScalar(Simulation* parent) :
	BondedPairArbitrary(parent) {
	m_function.setReturnType(Variant::SCALAR);
	m_datatype = DataFormat::DOUBLE;
	init();
}

BondedPairScalar::~BondedPairScalar() {
}

void BondedPairScalar::init() {
  // this sets name and className of m_properties
  m_properties.setClassName("ValCalculatorPair");
  // IMPORTANT: this sets back name
  m_properties.setName("BondedPairScalar");

  m_properties.setDescription("Module for calculation of a scalar Symbol stored per bonded pair of particles belonging to the specified list of bonded pairs.");

  m_expression = "1";  
  
}

//---- Methods ----

void BondedPairScalar::setup() {
  
  BondedPairArbitrary::setup();
 
}

