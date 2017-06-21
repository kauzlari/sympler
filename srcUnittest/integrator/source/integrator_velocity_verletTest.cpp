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



#include "integrator_velocity_verletTest.h"
#include <vector>


CPPUNIT_TEST_SUITE_REGISTRATION (IntegratorVelocityVerletTest);

void IntegratorVelocityVerletTest :: setUp (void)
{
  /*!
   * Initialize objects
   */
  particleTest = new Particle();
  simulationTest = new Simulation();
  controllerTest = new Controller(simulationTest);
  integratorVelocityVerletTest = new IntegratorVelocityVerlet(controllerTest);
}

void IntegratorVelocityVerletTest :: tearDown (void) 
{
  /*!
   * Delete objects
   */
  delete particleTest;
  delete simulationTest;
  delete controllerTest;
  delete integratorVelocityVerletTest;
}

void IntegratorVelocityVerletTest :: integrateVelocityTest (void)
{
  /*!
   * Test integrating velocity of a particle
   */
  particleTest->v = { { {0.0, 0.0, 0.0} } };
  particleTest->dt = 0.1;
  particleTest->force[0] = { {70.0, -2.0, 0.0} };
  integratorVelocityVerletTest->integrateVelocity(particleTest);
  point_t reference = { { {3.5, -0.1, 0.0} } };
  CPPUNIT_ASSERT_EQUAL (reference, particleTest->v);

  particleTest->v = { { {0.0, 0.0, 0.0} } };
  particleTest->dt = 0.01;
  particleTest->force[0] = { {-0.2, 2443, 90.0} };
  integratorVelocityVerletTest->integrateVelocity(particleTest);
  reference = { { {-0.001, 12.215, 0.45} } };
  CPPUNIT_ASSERT_EQUAL (reference, particleTest->v);

  particleTest->v = { { {0.0, 0.0, 0.0} } };
  particleTest->dt = 0.001;
  particleTest->force[0] = { {0.0, -87.34, 1.0} };
  integratorVelocityVerletTest->integrateVelocity(particleTest);
  reference = { { {0.0, -0.04367, 0.0005} } };
  CPPUNIT_ASSERT_EQUAL (reference, particleTest->v);
}

void IntegratorVelocityVerletTest :: hitPosTest (void)
{
  /*!
   * Test checking which of the times is the actual hit position
   */
  point_t hit_posTest = { { {0.0, 0.0, 0.0} } };
  particleTest->r = { { {-1.0, 2.0, 4.0} } };
  particleTest->v = { { {2.0, -2.0, 3.0} } };
  const point_t forceTest1 = { { {2.0, 2.0, 30.0} } };
  integratorVelocityVerletTest->hitPos(0.1, particleTest, hit_posTest, forceTest1);
  point_t reference = { { {-0.79, 1.81, 4.45} } };
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);

  hit_posTest = { { {0.0, 0.0, 0.0} } };
  particleTest->r = { { {0.0, 320.0, -1.4} } };
  particleTest->v = { { {-22.2, 5.1, 0.0} } };
  const point_t forceTest2 = { { {1.0, 0.0, 130.0} } };
  integratorVelocityVerletTest->hitPos(0.01, particleTest, hit_posTest, forceTest2);
  reference = { { {-0.22195, 320.051, -1.3935} } };
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);

  hit_posTest = { { {0.0, 0.0, 0.0} } };
  particleTest->r = { { {33.9, 20.0, -4.0} } };
  particleTest->v = { { {-1.2, 11.0, 0.0} } };
  const point_t forceTest3 = { { {0.0, 0.1, 0.55} } };
  integratorVelocityVerletTest->hitPos(0.001, particleTest, hit_posTest, forceTest3);
  reference = { { {33.8988, 20.01100005, -3.999999725} } };
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);
}


