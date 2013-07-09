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



#include <algorithm>

#include "gen_f.h"

#include "simulation.h"

using namespace std;

REGISTER_SMART_ENUM
(GenFType,
 "Force objects represent physical forces in the simulation. If more than one "
 "force is given the superposition principle is applied and the forces are "
 "summed. A force belongs to the Phase. A force always acts on exactly one colour " 
 "pair. This one has to be defined by the attributes 'species1' and 'species2'. "
 "Exception: One-particle forces (FParticleVels, FCentrifugal, FCoriolis, etc.).");


//---- Constructors/Destructor ----

GenF::GenF()
{
	throw gError("GenF::GenF(default)", "Should not be called. Contact the programmer.");
}


GenF::GenF(Simulation *simulation): Node((Node*) simulation)
{
	init();
}


GenF::~GenF()
{
}


//---- Functions ----

void GenF::init()
{
    STRINGPC
    (forceName, m_force_name,
     "Defines the force name to be used.");
     
  m_force_name = "default";
#ifdef _OPENMP
  m_offsetToVec.resize(global::n_threads);
//   m_posInVec = 0;
#endif
}

group_t group_intersection(const group_t &gr1, const group_t &gr2)
{
    if (gr1.empty())
        return gr2;
    else if (gr2.empty())
        return gr1;
    else {
        group_t intersection;

        set_intersection(
            gr1.begin(), gr1.end(),
            gr2.begin(), gr2.end(),
            inserter(intersection, intersection.begin())
        );

        return intersection;
    }
}
