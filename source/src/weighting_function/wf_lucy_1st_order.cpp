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



#include "wf_lucy_1st_order.h"


/* Register these WeightingFunctions with the factory. */
const WeightingFunction_Register<LucyWithWall1stOrder> lucy_with_wall_1st_order("LucyWithWall1stOrder");


/* --- LucyWithWall1stOrder --- */

/* Constructor/Destructor */

LucyWithWall1stOrder::LucyWithWall1stOrder(Node *parent): WeightingFunctionWithWall(parent) {
  init();
}


LucyWithWall1stOrder::~LucyWithWall1stOrder()
{
}


/* Methods */

void LucyWithWall1stOrder::init() {
  m_properties.setClassName("LucyWithWall1stOrder");
  m_properties.setDescription(
      "The Lucy interpolation function with 1st order correction at boundaries");
}


void LucyWithWall1stOrder::setup() {
  m_rc7 = pow(m_cutoff, 7);
  m_factor = 105/(16*M_PI*m_rc7);

  MSG_DEBUG("LucyWithWall1stOrder::setup", "m_rc7 = " << m_rc7);
}


