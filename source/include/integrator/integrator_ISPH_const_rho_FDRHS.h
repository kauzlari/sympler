/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


#ifndef __INTEGRATOR_ISPH_CONST_RHO_FDRHS_H
#define __INTEGRATOR_ISPH_CONST_RHO_FDRHS_H

#include "integrator_ISPH_const_rho.h"

using namespace std;

//----IntegratorIISPHconstRhoFDRHS ----

/*!
 * Child of abstract class \a IntegratorIISPHconstRho using the RHS
 * RHS = (rho_0-rho_adv)/rho_0/dt^2
 * where rho_0 is the constant reference density, rho_adv is the
 * advected density resulting from all forces except the pressure
 * gradient forces, and dt is the integration time step.
 */
class IntegratorISPHconstRhoFDRHS: public IntegratorISPHconstRho
{

 protected:
  
  /*!
   * Initialize the property list
   */
  void init();
  
  /*!
   * Adds the RHS of the PPE to the storage of the new pressure in the
   * iteration. In this class, 
   * RHS = (\a m_rho0 - rho_adv) / (\a m_rho0 * dt^2)
   * @param colour Colour of the particles this function should operate 
   * on
   */
  virtual void addRHStoNewPressure(size_t colour);


public:

  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorISPHconstRhoFDRHS(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorISPHconstRhoFDRHS();

  
};

#endif
