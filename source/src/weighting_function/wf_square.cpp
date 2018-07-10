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



#include "wf_square.h"


/* Register these WeightingFunctions with the factory. */
const WeightingFunction_Register<Square> square("Square");

/* --- Square --- */

/* Constructor/Destructor */

Square::Square(Node *parent): WeightingFunction(parent) {
  init();
}


Square::~Square()
{
}


/* Methods */

void Square::init() {
  m_properties.setClassName("Square");
  m_properties.setDescription(
      "Quadratic interpolation function:\nW(r)=(15/(2*pi*rc^3))*(1-r/rc)^2, with the cutoff distance 'rc'");
}


void Square::setup() {
  m_factor = 15/(2*M_PI*pow(m_cutoff, 5)); 
}

