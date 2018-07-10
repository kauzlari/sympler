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



#include "f_rand_var_noise.h"

#include "simulation.h"

#include "val_calculator_r_i.h"


using namespace std;

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

FRandVarNoise::FRandVarNoise(Simulation *simulation, double co, double friction,
  pair<IntegratorEnergy*, IntegratorEnergy*> ie, pair<string, string> species)
  : GenFSolo(simulation)
{
  init();

  m_factor = 4*friction;
  m_r_sqrt_dt = 1/sqrt(simulation->controller()->dt());
    
  m_cutoff = co;
  m_rcinv = 1/m_cutoff;

  m_ie = ie;

  m_species = species;

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

//   m_offset = m_cp->tagFormat().addAttribute
//     (FORCE_FACTOR_STR + className(), DataFormat::DOUBLE, false).offset;

  m_cp->setNeedPairs(true);
	
  m_cp->registerCalc(m_compute_ri_offset, new ValCalculatorRi, false);
  m_cp->setCutoff(m_cutoff);

  m_noise_offset = m_cp->tagFormat().addAttribute("noise", DataFormat::DOUBLE, false).offset;
}


FRandVarNoise::~FRandVarNoise()
{
}



//---- Methods ----

void FRandVarNoise::init()
{
  m_properties.setClassName("FRandVarNoise");

  m_properties.setDescription("Force using white noise as the proportional factor.");

  GenFSolo::init();
}



void FRandVarNoise::setup()
{
  GenFSolo::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_noise_offset = m_cp->tagFormat().addAttribute("noise", DataFormat::DOUBLE, false).offset;
}
