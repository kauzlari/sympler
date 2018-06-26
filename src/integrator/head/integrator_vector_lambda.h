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


#ifndef __INTEGRATOR_VECTOR_LAMBDA_H
#define __INTEGRATOR_VECTOR_LAMBDA_H

#include <math.h>

#include "integrator_vector.h"

using namespace std;

/*!
 * Adds an additional vector degree of freedom,
 * to the particles specified and integrates it with a second
 * order accurate predictor-corrector scheme.
 * Predictor step: predicted = old + lambda*dt*f_new
 * Corrector step: new = predicted + 0.5*dt*(f_new - f_old)
 * Usually, the forces are updated between the two steps.
 */
class IntegratorVectorLambda: public IntegratorVector
{
 protected:
  
  /*!
   * Factor for choice of predictor-corrector scheme.
   */
  double m_lambda;
  
  /*!
   * Helper being set to 0.5 - \a m_lambda in \a setup()
   */
  double m_lambda_diff;

  /*!
   * Initialize the property list
   */
  void init();
  
 public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */  
  IntegratorVectorLambda(Controller *controller);
  
  /*!
   * Destructor
   */
  virtual ~IntegratorVectorLambda();
  
  /*!
   * Register the field and the force of the field with the \a Particle
   */
  virtual void setup();
  
  /*!
   * Integrate the field
   */
  virtual void integrateStep1();
  
  /*!
   * Does nothing, because the predictor and corrector steps are merged
   */
  virtual void integrateStep2();

};

#endif
