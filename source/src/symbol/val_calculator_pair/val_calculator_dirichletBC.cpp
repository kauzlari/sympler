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


#include "val_calculator_dirichletBC.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"

// #include <utility>

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorDirichletBC::ValCalculatorDirichletBC(string symbol)
  : ValCalculatorBC(symbol)
{
//   MSG_DEBUG("ValCalculatorDirichletBC::ValCalculatorDirichletBC", "CONSTRUCTOR");
}

ValCalculatorDirichletBC::ValCalculatorDirichletBC(/*Node*/Simulation* parent)
  : ValCalculatorBC(parent)
{
  // next already done in upper class
  //  m_stage = 0;
  init();
}

void ValCalculatorDirichletBC::init()
{
  m_properties.setClassName("DirichletBCVels");

  m_wallSpecies = "undefined";
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorDirichletBC::setup()
{
//   M_CONTROLLER->registerForSetupAfterParticleCreation(this);
//   if(m_species.first == "undefined")
//     throw gError("ValCalculatorDirichletBC::setup", "Attribute 'species1' has value \"undefined\" .");
//   if(m_species.second == "undefined")
//     throw gError("ValCalculatorDirichletBC::setup", "Attribute 'species2' has value \"undefined\" .");
//   ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

//   // check whether wall particles are first or second in the pairs
//   if(m_wallSpecies == cp->firstSpecies())
//     m_wallIsSecond = false;
//   else if(m_wallSpecies == cp->secondSpecies())
//     m_wallIsSecond = true;
//   else
//     throw gError("ValCalculatorDirichletBC::setup", "Attribute wallSpecies = \"" + m_wallSpecies + "\" is neither equal to attribute species1 = \"" + m_species.first + "\" nor to species2 = \"" + m_species.second + "\".");

//   m_wallColour = M_MANAGER->getColour(m_wallSpecies);

  // m_datatype = DataFormat::POINT;

  ValCalculatorBC::setup();

}

// void /*pair<size_t, size_t>*/ ValCalculatorDirichletBC::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
// {
//   MSG_DEBUG("ValCalculatorDirichletBC::setSlot", "CALLED");
//   m_slot = slot = cp->tagFormat().addAttribute
//       ("ValCalculator_" + myName() + "_" + cp->toString(), DataFormat::POINT, false, m_symbolName).offset;
// }


// Commented out because this is just the specific velocity case, hence useless in an abstract parent class
#if 0

#ifndef _OPENMP
void ValCalculatorDirichletBC::compute(Pairdist* pD)
#else
void ValCalculatorDirichletBC::compute(Pairdist* pD, int thread_no)
#endif
{
  double innerDist;
  double outerDist;

  computeDists(pD, innerDist, outerDist);

    Particle* p1st = pD->firstPart();
    Particle* p2nd = pD->secondPart();

// //          point_t hitpos;
//     innerDist = HUGE_VAL;

     if(m_wallIsSecond)
     {

     }
     else // so !m_wallIsSecond
     {


          // if all worked fine, we may now compute the value of the outer particle
          // the check for innerDist != HUGE_VAL is done below
          // the first term is for Dirichlet != 0. We assume the value was assigned
          // to the wall particles.
        pD->tag.pointByOffset(m_slot) =
            (outerDist+innerDist)*(p1st->v)/innerDist-(outerDist/innerDist)*p2nd->v;
     } // of else of if(m_wallIsSecond)

    if(innerDist == HUGE_VAL)
      throw gError("ValCalculatorDirichletBC::compute", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + "c1=" + ObjToString(pD->firstPart()->c) + ", c2=" + ObjToString(pD->secondPart()->c) + ", r1=" + ObjToString(pD->firstPart()->r) + ", r2=" + ObjToString(pD->secondPart()->r));
         
//     //   } // this is the closing of the commented out if(m_cutoff > ...)                                                                      
}

#endif // end of #if 0
