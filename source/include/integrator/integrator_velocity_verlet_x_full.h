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


#ifndef __INTEGRATOR_VELOCITY_VERLET_X_FULL_H
#define __INTEGRATOR_VELOCITY_VERLET_X_FULL_H

#include "integrator_velocity_verlet_x.h"

using namespace std;

//----IntegratorVelocityVerletX ----

/*!
 * Modified Velocity-Verlet integrator for the positions and velocities
 * See: R. D. Groot and P. B. Warren, J. Chem. Phys. 107, 4423-4435 (1997).
 * For the integration of the positions, this Integrator uses a 
 * user-defined velocity. This velocity is also the one, which is integrated.
 */

class IntegratorVelocityVerletXFull: public IntegratorVelocityVerletX
{
  protected:
  
  public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
    IntegratorVelocityVerletXFull(Controller *controller);

  /*!
     * Destructor
   */
    virtual ~IntegratorVelocityVerletXFull();

    /*!
     * Initialise the property list
     */
    virtual void init();
    
    /*!
     * Setup this Integrator! Notice that this Integrator adds a velocity 
     * attribute same as \a IntegratorVelocityVerletX, BUT differently 
     * to \a IntegratorVelocityVerletX we make it persistent, because 
     * it is integrated !!!
     */
    virtual void setup();
    
    /*!
     * Correction of the velocity
     */
    virtual void integrateStep2();

    /*!
     * Prediction of the velocity
     */
    virtual void integrateVelocity(Particle* p);

};


#endif
  

