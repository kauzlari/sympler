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



#include "wf_input.h"


/* Register these WeightingFunctions with the factory. */
const WeightingFunction_Register<InputWF> InputWF("InputWF");


/* --- InputWF --- */

/* Constructor/Destructor */

InputWF::InputWF(Node *parent): WeightingFunction(parent) {
  init();
}


/*
InputWF::InputWF(double rc): WeightingFunction(rc) {
  init();
}
*/


InputWF::~InputWF()
{
}


/* Methods */

void InputWF::init() {
  m_properties.setClassName("InputWF");

  m_properties.setDescription(
      "Interpolant fully specified by the user. The whole interpolation function is given by W(r) = (1/Iw)*w(r), where w(r) is the shape of the interpolation function with w(0) = 1, and Iw is the area or volume integral of w(r) over its range. w(r) is specified by the attribute 'interpolation while Iw is specified by the attribute 'selfContribution'.");
    
  FUNCTIONFIXEDPC
      (interpolation , m_weightFunct, "The shape w(r) of the interpolant with w(0) = 1. The self contribution for r = 0 must be specified in 'selfContribution'."
  "'r' is the only accepted variable and represents the scalar distance between the particles");
  m_weightFunct.addVariable("r");
  //  m_scalar_name = "radius";
  //  m_scalar_symbol = "r";
  
  
  FUNCTIONFIXEDPC
    (selfContribution , m_selfContrib, "Prefactor 'Iw' of the weighting function. "
  "It is used for calculating self contributional values, i.e., for r = 0. "
      "'rc' is the only accepted variable and is a placeholder for the value defined for 'cutoff'.");
  m_selfContrib.addVariable("rc");
  
/*  FUNCTIONFIXEDPC
    (derivativeSC , m_Deriv_sC, "Prefactor of the derivative of the  weighting function"
  "This one is used for calculating self contributional values if r = NULL "
  "The variable rc is the cutoff radius");
  m_Deriv_sC.addVariable("rc");*/
  
  FUNCTIONFIXEDPC
      (weight , m_Deriv_wF, "This is -W'(r)/r, where W'(r) is the derivative of W(r) as defined above. I.e., the derivative has to be taken for the product '(1/Iw)*w(r)'. The user himself is responsible for consistency."
      "'r' is the only accepted variable and represents the scalar distance between the particles.");
  m_Deriv_wF.addVariable("r");
  
}


void InputWF::setup() {
 
}


