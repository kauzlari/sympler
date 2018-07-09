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



#ifndef __BOUNDARY_STEP_STOCH_H
#define __BOUNDARY_STEP_STOCH_H

#include "cell.h"
#include "boundary_arbitrary.h"

/*!
 * Simple step boundary with an inlet and an outlet
 */
class BoundaryStep: public BoundaryArbitrary
{
protected:
  double m_big_height;
  double m_big_lflow;
  double m_before_step;
  double m_after_step;
  double m_small_lflow, m_small_height, m_width;
  int m_flow_dir, m_step_loc, m_third_dir;
  bool m_side_periodic;
  string m_side_wall;
  double m_temperature_sides, m_temperature_step_dir;
    
  void init();
	
public:
  BoundaryStep(Phase *phase);
  virtual ~BoundaryStep();
	
  virtual void setup(Simulation* sim, ManagerCell *mgr);

  virtual void setup();
};

#endif
