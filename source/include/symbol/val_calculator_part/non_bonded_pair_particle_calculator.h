/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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



#ifndef __NON_BONDED_PAIR_PARTICLE_CALCULATOR_H
#define __NON_BONDED_PAIR_PARTICLE_CALCULATOR_H

#include "symbol.h"
#include "val_calculator_pair.h"
#include "val_calculator_part.h"

// #define VC_MAX_STAGE 1

class Pairdist;
class ColourPair;


class NonBondedPairParticleCalculator : public ValCalculatorPart
{

 protected:

/*!
 * Should all colour combinations be considered. This disables \a m_species.
 */
  bool m_allPairs;


  /*!
   * initialise this Caluclator
  */
  virtual void init();

  public:

  /*!
   * Constructor
   */
  NonBondedPairParticleCalculator(string symbol);


  /*!
     * Constructor for Node hierarchy
   */
    NonBondedPairParticleCalculator(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~NonBondedPairParticleCalculator() {
    }

};


#endif
