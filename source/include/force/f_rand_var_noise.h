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



#ifndef __F_RAND_VAR_NOISE_H
#define __F_RAND_VAR_NOISE_H 

#include "random.h"
#include "gen_f_solo.h"
#include "integrator_energy.h"

using namespace std;


//---- Frand ----

/*!
 * Random force with variable noise amplitude. Fixme!!! Probably obsolete.
 */
class FRandVarNoise : public GenFSolo
{
protected:
  double m_factor, m_r_sqrt_dt;
  pair<IntegratorEnergy*, IntegratorEnergy*> m_ie;

  int m_noise_offset;
    
  void init();

  inline virtual double computeForceFactor(Pairdist *pair) {
  
   // RandomNumberGenerator m_rng;
  
    double noise = sqrt(m_factor / 
      (m_ie.first->reciprocalTemperature(*pair->firstPart()) +
       m_ie.second->reciprocalTemperature(*pair->secondPart())));

    pair->tag.doubleByOffset(m_noise_offset) = noise;

    return
      m_rng.normal(1) * noise * m_r_sqrt_dt
      * (pair->tag.doubleByOffset(m_compute_ri_offset) - m_rcinv);
  }

public:
  FRandVarNoise(Simulation *simulation, double co, double friction,
                pair<IntegratorEnergy*, IntegratorEnergy*> ie, pair<string, string> species);
  virtual ~FRandVarNoise();

  virtual void setup();
};

#endif
