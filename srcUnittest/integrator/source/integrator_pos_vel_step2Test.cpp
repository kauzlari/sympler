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


#include "integrator_pos_vel_step2Test.h"
#include <vector>


CPPUNIT_TEST_SUITE_REGISTRATION (IntegratorPosVelStep2Test);

void IntegratorPosVelStep2Test :: setUp (void)
{
  IntegratorPositionTest :: setUp();
  
  m_integratorPositionTest = new IntegratorPosVelStep2(m_controllerTest);
}

void IntegratorPosVelStep2Test :: tearDown (void) 
{
  IntegratorPositionTest :: tearDown();
  
  delete m_integratorPositionTest;
}

void IntegratorPosVelStep2Test :: integrateVelocityTest (void)
{
  IntegratorPositionTest :: integrateVelocityTest();
}

void IntegratorPosVelStep2Test :: hitPosTest (void)
{
  point_t hit_posTest = { { {0.0, 0.0, 0.0} } };
  m_particleTest->r = { { {-1.0, 2.0, 4.0} } };
  m_particleTest->v = { { {2.0, -2.0, 3.0} } };
  point_t forceTest = { { {2.0, 2.0, 30.0} } };
  m_integratorPositionTest->hitPos(0.1, m_particleTest, hit_posTest, forceTest);
  point_t reference = { { {-0.8, 1.8, 4.3} } };
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);

  hit_posTest = { { {0.0, 0.0, 0.0} } };
  m_particleTest->r = { { {0.0, 320.0, -1.4} } };
  m_particleTest->v = { { {-22.2, 5.1, 0.0} } };
  forceTest = { { {1.0, 0.0, 130.0} } };
  m_integratorPositionTest->hitPos(0.01, m_particleTest, hit_posTest, forceTest);
  // yes, same z-value because force does not have an effect and vz = 0!
  reference = { { {-0.222, 320.051, -1.4} } };
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);

  hit_posTest = { { {0.0, 0.0, 0.0} } };
  m_particleTest->r = { { {33.9, 20.0, -4.0} } };
  m_particleTest->v = { { {-1.2, 11.0, 0.0} } };
  forceTest = { { {0.0, 0.1, 0.55} } };
  m_integratorPositionTest->hitPos(0.001, m_particleTest, hit_posTest, forceTest);
  reference = { { {33.8988, 20.011, -4.0} } };
  CPPUNIT_ASSERT_EQUAL (reference, hit_posTest);
}


