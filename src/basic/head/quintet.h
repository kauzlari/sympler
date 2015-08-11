/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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


#ifndef __QUINTET_H
#define __QUINTET_H

#include "misc.h"
#include "quintet_calculator.h"


/*!
 * Used for storing bonded quintets of \a Particle s in the \a Phase
 */
struct quintet_t{
  /*!
   * Particles defined from r00
   */
  /*!
   * first particle
   */
  Particle* p00;
  /*!
   * second particle
   */
  Particle* p20;
  /*!
   * third particle
   */
  Particle* p22;
  /*!
   * fourth particle
   */
  Particle* p02;
  /*!
   * center particle
   */
  Particle* p11;

  
  /*!
   * The \a Phase, this \a quintet_t belongs to. Currently (2015-07-20) there is only one \a Phase in the \a Simulation
   */
  static Phase* s_phase;
  
  /*!
   * As the name says
   */
  void runBondedQuintetCalculators(size_t stage, size_t listIndex)
  {
    FOR_EACH
      (vector<QuintetCalculator*>,
       s_phase->bondedQuintetCalculators(stage, listIndex),
       (*__iFE)->compute(this);
       );
  }
  
  /*!
   * As the name says
   */
  void runBondedQuintetCalculators_0(size_t stage, size_t listIndex)
  {
    FOR_EACH
      (vector<QuintetCalculator*>,
       s_phase->bondedQuintetCalculators_0(stage, listIndex),
       (*__iFE)->compute(this);
       );
  }
  
};

#endif
