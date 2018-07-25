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


#ifndef __INTEGRATOR_SCALAR_STEP2_H
#define __INTEGRATOR_SCALAR_STEP2_H


#include <math.h>

#include "integrator_scalar.h"

using namespace std;


/*!
 * Adds an additional scalar degree of freedom, to the particles of 
 * the specified colour and integrates it according to the following 
 * scheme:
 * integration-step1: no activity 
 * integration-step2: s(t + dt) = s(t) + dt * F(t)
 * Here, s is the scalar, F is its flux (usually computed with Force 
 * modules such as FPairScalar or FParticleScalar), t is time and dt is 
 * the size of the integration time step (defined in the Controller). 
 * Further information on the integration-steps, including their place 
 * in the total SYMPLER workflow, can be found with the help option 
 * "--help workflow".
 */
class IntegratorScalarStep2: public IntegratorScalar
{
  protected:
  
  /*!
   * Factor for choice of predictor-corrector scheme.
   * FIXME: Every Integrator*Lambda currently has this one 
   * -> generalise in a parent
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
  IntegratorScalarStep2(Controller *controller);
  
  /*!
   * Destructor
   */
  virtual ~IntegratorScalarStep2();
  
  /*!
   * Register the field and the force of the field with the \a Particle
   */
  virtual void setup();
  
  /*!
   * Does nothing in the scheme implemented here
   */
  virtual void integrateStep1();
  
  /*!
   * Integration s(t + dt) = s(t) + dt * F(t)
   */
  virtual void integrateStep2();

};

#endif
