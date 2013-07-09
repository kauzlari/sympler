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



#include "f_rand.h"
#include "simulation.h"
// #include "colour_pair.h"

#include "val_calculator_r_i.h"

using namespace std;

const GenFTypeConcr<Frand> frand("Frand");

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

Frand::Frand(Simulation *simulation)
    : GenFSolo(simulation)
{
    init();
}


Frand::Frand(Simulation *simulation, double co, double ns, pair<string, string> species)
  : GenFSolo(simulation)
{
  init();
    
  m_noise = ns;
  m_noise_and_time = m_noise/sqrt(simulation->controller()->dt());
    
  MSG_DEBUG("Frand::Frand", "noise = " << m_noise);

  m_cutoff = co;
  m_rcinv = 1/m_cutoff;
  m_species = species;

  ColourPair *cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  
//   m_offset = m_cp->tagFormat().addAttribute
//     (FORCE_FACTOR_STR + className(), DataFormat::DOUBLE, false).offset;

  cp->setNeedPairs(true);
	
  cp->registerCalc(m_compute_ri_offset, new ValCalculatorRi, false);
  cp->setCutoff(m_cutoff);
  
//   m_interpolant.setColourPair(m_cp);
}


Frand::~Frand()
{
}



//---- Methods ----

void Frand::init()
{
  m_properties.setClassName("Frand");

  m_properties.setDescription("Force implementing white pair-wise noise. Note: This force already uses a linear weighting function of the form (1-rij/rc) where rij is the particle distance and rc the cutoff. If you require a normalisation factor, then include it in the 'noise' attribute. If you wish to use a different weighting function, you can modify it by using the WeightingFuncton module InputWF with non-trivial entries. Frand multiplies with the 'interpolate'-entry of the used WeightingFunction.");
    
  DOUBLEPC
    (noise, m_noise, 0,
     "Noise amplitude. Should be >0");
    
  m_noise = 1;

//     GenFSolo::init();
}

/*virtual*/ void Frand::setup()
{
	GenFSolo::setup();

	m_noise_and_time = m_noise/sqrt(M_SIMULATION->controller()->dt());
}

