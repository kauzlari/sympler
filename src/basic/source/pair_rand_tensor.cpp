/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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


#include "pair_rand_tensor.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"

using namespace std;

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

const SymbolRegister<PairRandTensor> pair_rand_tensor("PairRandTensor");

PairRandTensor::PairRandTensor(string symbol) :
	PairRandArbitrary(symbol) {
}

PairRandTensor::PairRandTensor(Simulation* parent) :
	PairRandArbitrary(parent) {
	m_function.setReturnType(Variant::TENSOR);
	m_datatype = DataFormat::TENSOR;
	init();
}

PairRandTensor::~PairRandTensor() {
}

void PairRandTensor::init() {

  m_properties.setClassName("PairRandTensor");
  m_properties.setDescription("Module for calculation of a random square matrix Symbol [pairFactor]*[ranTensor] stored per non-bonded pair of particles, where [pairFactor] is a user-defined runtime-compiled suqare matrix expression and [ranTensor] is a square matrix of independent Gaussian random numbers with zero mean and unit variance. Note that you can modify the mean and variance through a suitable [pairFactor]. [pairFactor] and [ranTensor] are multiplied componentwise.");
  
  m_expression_string = "unitMat(1)";

}
