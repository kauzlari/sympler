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



#include "manager_cell.h"
#include "simulation.h"
#include "weighting_function.h"
#include "colour_pair.h"

#include "f_pair_arbitrary_wf.h"


#define M_SIMULATION ((Simulation *) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

FPairArbitraryWF::FPairArbitraryWF(Simulation *simulation): FPairArbitrary(simulation)
{
  init();
}


FPairArbitraryWF::~FPairArbitraryWF()
{
}


void FPairArbitraryWF::init()
{
  m_properties.setName("FPairArbitraryWF");

  m_properties.setDescription("Pair force with associated weighting functon.");
    


  STRINGPC
    (weightingFunction, m_weighting_function,
     "Defines the weighting function to be used. Note that each expression, which is used"
     " in this force and needing a cutoff, and for which no explicit cutoff is specified, "
     "uses the cutoff defined in the chosen weighting function.");
     
  m_weighting_function = "default";

}


void FPairArbitraryWF::setup()
{
  FPairArbitrary::setup();
  
  m_wf = M_SIMULATION->findWeightingFunction(m_weighting_function);
  m_cutoff = m_wf->cutoff();

  ColourPair *cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  	
  cp->setCutoff(m_cutoff);
}

