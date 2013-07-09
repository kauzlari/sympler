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



#ifndef __F_RAND_H
#define __F_RAND_H

#include "random.h"
#include "gen_f_solo.h"

using namespace std;


//---- Frand ----

/*!
 * Random force
 */
class Frand : public GenFSolo
{
protected:
  /*!
   * The noise amplitude
   */
  double m_noise;

  /*!
   * = m_noise/sqrt(dt)
   */
  double m_noise_and_time;
    
  void init();

  inline virtual double computeForceFactor(Pairdist *pair) {
     return m_rng.normal(1) * m_noise_and_time * (1 - m_rcinv*pair->abs());
  }
                                             
public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  Frand(Simulation *simulation);

  /*!
   * Constructor. Used from FFlukDiss to construct the fluctuating force.
   * Fixme!!! Probably obsolete!
   * @param simulation Pointer to the simulation object
   * @param co Cut-off radius
   * @param ns Noise amplitude
   * @param species Species, this force should act on
   */
  Frand(Simulation *simulation, double co, double ns, pair<string, string> species);
  virtual ~Frand();

  /*!
   * Return the noise amplitude
   */
  double noise() const {
    return m_noise;
  }

  virtual void setup();
};

#endif
