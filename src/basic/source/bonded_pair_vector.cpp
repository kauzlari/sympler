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


#include "bonded_pair_vector.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

const SymbolRegister<BondedPairVector> bonded_pair_vector("BondedPairVector");

BondedPairVector::BondedPairVector(Simulation* parent) :
	BondedPairArbitrary(parent) {
	m_function.setReturnType(Variant::VECTOR);
	m_datatype = DataFormat::POINT;
	init();
}

BondedPairVector::~BondedPairVector() {
}

void BondedPairVector::init() {
  // this set name and className of m_properties
  m_properties.setClassName("ValCalculatorPair");
  // IMPORTANT: this sets back name for correct identification of plymorphic type
  m_properties.setName("BondedPairVector");

  m_properties.setDescription("Module for calculation of a vectorial Symbol stored per bonded pair of particles belonging to the specified list of bonded pairs.");

  m_expression = "idVec(1)";  
  
}

//---- Methods ----

void BondedPairVector::setup() {
  
  BondedPairArbitrary::setup();
 
}

// Now (2014-10-31) in BondedPairArbitrary
// void BondedPairVector::setSlot(ColourPair* cp, size_t& slot,
// 						 bool oneProp) {
//   m_slot = slot = cp->tagFormat().addAttribute
//     // OLD-STYLE. No idea why this was done.
//     //    ("BondedPairVector_" + cp->toString(), DataFormat::DOUBLE, false, "m_scalar").offset;
//     // NEW-STYLE (2014-10-31)
//     (m_symbolName, DataFormat::POINT, false, m_symbolName).offset;
// }

