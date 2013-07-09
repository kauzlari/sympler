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


#ifndef __TRIPLET_H
#define __TRIPLET_H 

#include "misc.h"
#include "triplet_calc_angular_f.h"

/*!
 * Used for storing bonded triplets of \a Particle s in the \a Phase
 */
struct triplet_t{
  /*!
   * first particle
   */
  Particle* a;
  /*!
   * second particle
   */
  Particle* b;
  /*!
   * third particle
   */
  Particle* c;
  
  /*!
   * The \a Phase, this \a triple_t belongs to. Currently (2009-07-03) there is only one \a Phase in the \a Simulation
   */
  static Phase* s_phase;
  
  /*!
   * As the name says
   */
  void runBondedTripletCalculators(size_t stage, size_t listIndex)
  {
    FOR_EACH
      (vector<TripletCalculator*>,
       s_phase->bondedTripletCalculators(stage, listIndex),
       (*__iFE)->compute(this);
       );
  }
  
  /*!
   * As the name says
   */
  void runBondedTripletCalculators_0(size_t stage, size_t listIndex)
  {
    FOR_EACH
      (vector<TripletCalculator*>,
       s_phase->bondedTripletCalculators_0(stage, listIndex),
       (*__iFE)->compute(this);
       );
  }
  
};

#endif
