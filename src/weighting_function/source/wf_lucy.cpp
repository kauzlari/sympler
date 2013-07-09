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



#include "wf_lucy.h"


/* Register these WeightingFunctions with the factory. */
const WeightingFunction_Register<Lucy> lucy("Lucy");


/* --- Lucy --- */

/* Constructor/Destructor */

Lucy::Lucy(Node *parent): WeightingFunction(parent) {
  init();
}


/*
Lucy::Lucy(double rc): WeightingFunction(rc) {
  init();
}
*/


Lucy::~Lucy()
{
}


/* Methods */

void Lucy::init() {
  m_properties.setClassName("Lucy");
  m_properties.setDescription(
      "The Lucy interpolation function:\nW(r)=(105/(16*pi*rc^3))*(1+3*r/rc)*(1-r/rc)^3, with the cutoff distance 'rc'");
}


void Lucy::setup() {
  m_factor_i = 105/(16*M_PI*pow(m_cutoff, 7));
  m_factor_w = 315/(4*M_PI*pow(m_cutoff, 7));
}


