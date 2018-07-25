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


#include "integrator_positionTest.h"
#include <vector>


void IntegratorPositionTest :: setUp (void)
{

  // Initialize objects
  m_particleTest = new Particle();
  m_simulationTest = new Simulation();
  m_controllerTest = new Controller(m_simulationTest);
}

void IntegratorPositionTest :: tearDown (void) 
{
  // Delete objects
  delete m_controllerTest;
  delete m_simulationTest;
  delete m_particleTest;
}

void IntegratorPositionTest :: integrateVelocityTest (void)
{
  // Here we test a function that should do nothing
  point_t reference = m_particleTest->v = { { {3.5, -0.1, 0.0} } };
  m_integratorPositionTest->integrateVelocity(m_particleTest);
  CPPUNIT_ASSERT_EQUAL (reference, m_particleTest->v);

}

void IntegratorPositionTest :: hitPosTest (void)
{
  // Here we test a function that should do nothing
  point_t reference = { { {-0.79, 1.81, 4.45} } };
  point_t hit_posTest = reference;
  m_particleTest->r = { { {-1.0, 2.0, 4.0} } };
  m_particleTest->v = { { {2.0, -2.0, 3.0} } };
  const point_t forceTest1 = { { {2.0, 2.0, 30.0} } };
  m_integratorPositionTest->hitPos(0.1, m_particleTest, hit_posTest, forceTest1);
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);
}


