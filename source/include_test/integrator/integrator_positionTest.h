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



#ifndef INTEGRATOR_POSITION_TEST_H
#define INTEGRATOR_POSITION_TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "integrator_position.h"
#include "particle.h"
#include "simulation.h"
#include "controller.h"


using namespace std;

class IntegratorPositionTest : public CPPUNIT_NS :: TestFixture
{

 public:
  
  /*!
   * Initialize objects
   */
  void setUp (void);

  /*!
   * Delete objects
   */
  void tearDown (void);

  
 protected:

  /*!
   * Test integrating velocity of a particle
   */
  void integrateVelocityTest (void);

  /*!
   * Test computation of hit position
   */
  void hitPosTest (void);

  /*!
   * Dummy particle for testing
   */
  Particle *m_particleTest;

  /*!
   * Instance of child of \a IntegratorPosition to be tested
   */
  IntegratorPosition *m_integratorPositionTest;

  /*!
   * Required to instantiate an \a IntegratorPosition
   */
  Controller *m_controllerTest;
  
 private:
  
  Simulation *m_simulationTest;
  
};

#endif
