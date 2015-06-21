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


#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "triplet_calc_part_scalar.h"

// #include "triplet.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
// #define PI 3.141592654
#define PI M_PI

TripletCalcPartScalar::TripletCalcPartScalar(Simulation *simulation): TripletCalculator(simulation)
{
  // For now (2015-06-18) all inheriting classes do not depend on other symbols
  // if this changes, setting m_stage must be moved down the hierarchy (or just reset it)
  m_stage = 0;

  m_datatype = DataFormat::DOUBLE;
  
  init();
#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void TripletCalcPartScalar::init()
{
  m_properties.setClassName("TripletCalcPart");
  m_properties.setName("TripletCalcPartScalar");

  FUNCTIONFIXEDPC
      (expression, m_expression, "Scalar function depending on the variable \"ca\" which is the cosine of the triplet-angle alpha.\nNOTE OF CAUTION: The following definitions are assumed: If i, j, k are, respectively, the \"left\", \"central\", and \"right\" particle of the triplet, then the vectors rij=ri-rj and rjk=rj-rk enclose an angle psi=pi-alpha for which cos(psi)=cos(pi-alpha)=-cos(alpha)=rij.rjk/(|rij||rjk|)=-\"ca\". ");

  m_expression.addVariable("ca");

}


void TripletCalcPartScalar::setup()
{
  TripletCalculator::setup();
}


void TripletCalcPartScalar::setupAfterParticleCreation()
{
  TripletCalculator::setupAfterParticleCreation();
}


/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool TripletCalcPartScalar::findStage()
{
  
  // currently (2010/05/17) this always returns true
  return TripletCalculator::findStage();
}

/*!
 * Determines \a m_stage of the current \a Symbol.
 * By default, we assume that the stage is fixed and known during compile-time, 
 * so this function does nothing except returning the message (true) that the 
 * stage was already found. Symbols, which determine the stage during run-time 
 * have to redefine this function.
 */
bool TripletCalcPartScalar::findStage_0()
{
  // currently (2010/05/17) this always returns true
  return TripletCalculator::findStage_0();
}
