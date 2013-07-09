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


#ifndef __INTEGRATOR_TENSOR_LAMBDA_H
#define __INTEGRATOR_TENSOR_LAMBDA_H


#include <math.h>

// #include "function.h"
// #include "particle.h"
#include "integrator_tensor.h"

using namespace std;

// class GenF;
// class Phase;
// class Controller;

/*!
 * Second order accurate integrator of a tensor field
 */
class IntegratorTensorLambda: public IntegratorTensor
{
  protected:
  
  /*!
   * Factor for choice of predictor-corrector scheme. Note that, currently,
   * the predictor and corrector steps are merged together, since it is 
   * assumed that only the predicted value is needed for computation of the 
   * other quantities
   */
    double m_lambda;

    /*! helper
     */
    double m_lambda_sum;
    
    /*! helper
     */
    double m_lambda_diff;
    
    /*!
     * Boolean helper for distinguishing between the first integration and 
     * later ones
     */
    bool m_laterStep;
    
  /*!
     * Initialize the property list
   */
    void init();

  public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */  
    IntegratorTensorLambda(Controller *controller);

  /*!
     * Destructor
   */
    virtual ~IntegratorTensorLambda();

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
