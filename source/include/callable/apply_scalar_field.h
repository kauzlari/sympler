/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#ifndef __APPLY_SCALAR_FIELD_H
#define __APPLY_SCALAR_FIELD_H

#include "thermostat.h"
#include "function_fixed.h"

using namespace std;

class Phase;
class Simulation;

/*!
 * \a Callable setting a user-defined scalar field to a user-defined 
 * scalar expression
 * FIXME: Put this class and \a ApplyVectorField (and \a ApplyVelField ? and 
 * others?) into a hierarchy with parent CblParticleArbitrary (DO NOT 
 * inherit from Thermostat anymore! Check \a CblPairParticleArbitrary .) and 
 * switch the \a FunctionFixed to the more flexible \a FunctionParticle. 
 */
class ApplyScalarField : public Thermostat
{

 protected:

  /*!
   * Sets the user-defined scalar specified by \a m_symbolName to the specified algebraic scalar expression. See \a init() for the known symbols.
   */
  FunctionFixed m_expr;
  
  /*!
   * The colour, this \a Callable should act on
   */
  size_t m_colour;
 
  /*!
   * The species, this \a Callable should act on
   */
  string m_species;
 
  /*!
   * The name of the \a Symbol, this \a Callable should act on
   */
  string m_symbolName;
 
  /*!
   * The offset of the \a Symbol, this \a Callable should act on
   */
  size_t m_offset;
 
  /*!
   * Initialize the property list
   */
  void init();

 public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  ApplyScalarField(Simulation* sim);

  /*!
   * Destructor
   */
  virtual ~ApplyScalarField() {}

  /*!
   * Thermalize the system
   */
  virtual void thermalize(Phase* p);

  /*!
   * Setup this \a Thermostat
   */
  virtual void setup();
};

#endif

