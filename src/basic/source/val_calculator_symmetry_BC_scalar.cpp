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


#include "val_calculator_symmetry_BC_scalar.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"
#include "particle_cache.h"

// #include <utility>

const SymbolRegister<ValCalculatorSymmetryBCScalar> val_calc_symmetry_BC_scalar("ValCalculatorSymmetryBCScalar");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorSymmetryBCScalar::ValCalculatorSymmetryBCScalar(string symbol)
  : ValCalculatorBC(symbol)
{
//   MSG_DEBUG("ValCalculatorSymmetryBCScalar::ValCalculatorSymmetryBCScalar", "CONSTRUCTOR");
}

ValCalculatorSymmetryBCScalar::ValCalculatorSymmetryBCScalar(/*Node*/Simulation* parent)
  : ValCalculatorBC(parent)
{
  // next already done in upper class; there is a findStage() for this class
  //  m_stage = 0;
  init();
}

void ValCalculatorSymmetryBCScalar::init()
{
  m_properties.setClassName("SymmetryBCScalar");

  m_properties.setDescription("Saves the pair-specific value of an arbitrary scalar of the boundary particle used for applying a symmetry boundary condition (BC) in each pair of particles. This calculator uses a linear approximation and relies on a local gradient of the scalar which has to be computed before. The gradient is assigned to this Calculator by specifying the attribute 'gradient'.");

  STRINGPC
      (scalar, m_scalarName,
       "Name of the scalar the symmetry BC should be computed for.");
  
  STRINGPC
      (gradient, m_gradientName,
       "Name of the gradient of the scalar to be used.");
  
  m_symbolName = "undefined";
  m_wallSpecies = "undefined";
  m_gradientName = "undefined";
  m_scalarName = "undefined";
  
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorSymmetryBCScalar::setup()
{
  if(m_symbolName == "undefined")
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Attribute 'symbol' has value \"undefined\" .");

  if(m_gradientName == "undefined")
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Attribute 'gradient' has value \"undefined\" .");

  if(m_scalarName == "undefined")
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Attribute 'scalar' has value \"undefined\" .");



  m_datatype = DataFormat::DOUBLE;
  
  ValCalculatorBC::setup();  

  // now, ValCalculatorBC::setup() has defined m_species
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  size_t freeColour;
  string freeSpecies;
  if(m_wallIsSecond)
    {
      freeColour = cp->firstColour();
      freeSpecies = cp->firstSpecies();
    }
  else
    {
      freeColour = cp->secondColour();
      freeSpecies = cp->secondSpecies();
    }

  if(!Particle::s_tag_format[freeColour].attrExists(m_gradientName))
        throw gError("ValCalculatorSymmetryBCScalar::setup", "Symbol " + m_gradientName + " not found for species '" + freeSpecies + "'.");    

  if(Particle::s_tag_format[cp->firstColour()].attrExists(m_symbolName))
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Symbol " + m_symbolName + " already exists for species '" + cp->firstSpecies() + "'.");
    
  if(Particle::s_tag_format[cp->secondColour()].attrExists(m_symbolName))
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Symbol " + m_symbolName + " already exists for species '" + cp->secondSpecies() + "'.");

  if(!Particle::s_tag_format[cp->firstColour()].attrExists(m_scalarName))
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Symbol " + m_scalarName + " not found for species '" + cp->firstSpecies() + "'.");
    
  if(!Particle::s_tag_format[cp->secondColour()].attrExists(m_scalarName))
    throw gError("ValCalculatorSymmetryBCScalar::setup", "Symbol " + m_scalarName + " not found for species '" + cp->secondSpecies() + "'.");
    
  m_gradientOffset = 
        Particle::s_tag_format[freeColour].offsetByName/*indexOf*/(m_gradientName);

  m_scalarOffset.first = 
    Particle::s_tag_format[cp->firstColour()].offsetByName/*indexOf*/(m_scalarName);

  m_scalarOffset.second = 
    Particle::s_tag_format[cp->secondColour()].offsetByName/*indexOf*/(m_scalarName);


}

// void /*pair<size_t, size_t>*/ ValCalculatorSymmetryBCScalar::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
// {
//   MSG_DEBUG("ValCalculatorSymmetryBCScalar::setSlot", "CALLED");
//   m_slot = slot = cp->tagFormat().addAttribute
//       ("ValCalculator_" + myName() + "_" + cp->toString(), DataFormat::POINT, false, m_symbolName).offset;
// }

// FIXME: inline ?
void ValCalculatorSymmetryBCScalar::compute(Pairdist* pD)
{
  double innerDist;
  double outerDist;

  computeDists(pD, innerDist, outerDist);

    Particle* p1st = pD->firstPart();
    Particle* p2nd = pD->secondPart();

    if(m_wallIsSecond)
      {
	// Compute the projection of the gradient onto the normal of the wall
	// The normal points inside, so the projected gradient always describes 
	// how the quantity changes from the wall to the inside.
	double gradProj = p1st->tag.pointByOffset(m_gradientOffset)*m_currentWall->normal();

        // if all worked fine, we may now compute the value of the outer particle
        // the check for innerDist != HUGE_VAL is done below
	pD->tag.doubleByOffset(m_slot) = 
	  p1st->tag.doubleByOffset(m_scalarOffset.first) + (outerDist-innerDist)*gradProj;


      }
    else
      {
	// Compute the projection of the gradient onto the normal of the wall
	// The normal points inside, so the projected gradient always describes 
	// how the quantity changes from the wall to the inside.
	double gradProj = p2nd->tag.pointByOffset(m_gradientOffset)*m_currentWall->normal();

        // if all worked fine, we may now compute the value of the outer particle
        // the check for innerDist != HUGE_VAL is done below
        // the first term is for Dirichlet != 0. We assume the value was assigned 
        // to the wall particles.
	pD->tag.doubleByOffset(m_slot) = 
	  p2nd->tag.doubleByOffset(m_scalarOffset.second) + (outerDist-innerDist)*gradProj;

//            MSG_DEBUG("ValCalculatorSymmetryBCScalar::compute", "WALLISFIRST: c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << ", tester=" << tester);
      }

    if(innerDist == HUGE_VAL)
      throw gError("ValCalculatorSymmetryBCScalar::compute", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + "c1=" + ObjToString(pD->firstPart()->c) + ", c2=" + ObjToString(pD->secondPart()->c) + ", r1=" + ObjToString(pD->firstPart()->r) + ", r2=" + ObjToString(pD->secondPart()->r));
         
//     //   } // this is the closing of the commented out if(m_cutoff > ...)             
}


// FIXME: inline ?
// FIXME: it is just because of the line m_currentWall = ... that we have to repeat the definition of this function !!! This line is needed for computing the projection of the gradient onto the normal of the wall
void ValCalculatorSymmetryBCScalar::computeDists(Pairdist* pD, double& innerDist, double& outerDist)
{
//   double innerDist;
//   double outerDist;
  double tempDist;
  vector<pair<Wall*, double> >* walls;
  vector<pair<Wall*, double> >::iterator wallIt;
  vector<pair<Wall*, double> >::iterator wallEnd;

  // it does not make sense to check for a cutoff because we should make sure that also the interaction with the longest range may safely use the boundary values computed here
/*  if (m_cutoff > pD->abs()) 
  {                                    */
//     point_t temp;
// 
//     point_t tester;
//     double arg[3];

    Particle* p1st = pD->firstPart();
    Particle* p2nd = pD->secondPart();
//          point_t hitpos;
    innerDist = HUGE_VAL;
    m_currentWall = NULL;
    if(m_wallIsSecond)
    {
           // find the wall-segment where the rij-vector is hitting and the distance of ri to it
      wallEnd = m_wallTable[p2nd->mySlot]->end();
      for(wallIt = m_wallTable[p2nd->mySlot]->begin(); wallIt != wallEnd; ++wallIt)
      {
        if(wallIt->first->intersects
           (
           /*this is either [ri] or if necessary its periodic image*/ 
           p2nd->r + pD->cartesian(), 
            /*const point_t&*/ -pD->cartesian(), // !!! -[rij] points from i to j !!!
                /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
           )
          )
        {
          if(tempDist < innerDist)
          {
            innerDist = tempDist; 
            outerDist = wallIt->second;
            m_currentWall = wallIt->first;
          }
        }
      }
      walls = m_periodicWallTable[p2nd->mySlot];
          // if it is NULL (was set in setupAfter... as initial value) we may skip it
      if(walls)
      {
        wallEnd = walls->end();
        for(wallIt = walls->begin(); wallIt != wallEnd; ++wallIt)
        {
          if(wallIt->first->intersects
             (
                /*here we may always take [ri] since we used the periodic image of the 
             wall particle for the computation of its distance to the wall*/ 
             p1st->r, 
                  /*const point_t&*/ -pD->cartesian(), // !!! -[rij] points from i to j !!!
                  /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
             )
            )
          {
            if(tempDist < innerDist)
            {
              innerDist = tempDist; 
              outerDist = wallIt->second;
	      m_currentWall = wallIt->first;
            }
          }
        }

      }




    }
    else // so !m_wallIsSecond
    {
           // find the wall-segment (-1)*rij-vector is hitting and the distance of rj to it
      wallEnd = m_wallTable[p1st->mySlot]->end();
      for(wallIt = m_wallTable[p1st->mySlot]->begin(); wallIt != wallEnd; ++wallIt)
      {
//               MSG_DEBUG("ValCalculatorSymmetryBCScalar::computeDists", "WALLISFIRST, non-periodic: checking, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
        if(wallIt->first->intersects
           (
           /*this is either [rj] or if necessary its periodic image*/ 
           p1st->r - pD->cartesian(), 
                /*const point_t&*/ pD->cartesian(), // !!! [rij] points from j to i !!!
                /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
           )
          )
        {
          if(tempDist < innerDist)
          {
            innerDist = tempDist; 
            outerDist = wallIt->second;
            m_currentWall = wallIt->first;
//                  MSG_DEBUG("ValCalculatorSymmetryBCScalar::computeDists", "WALLISFIRST, non-periodic: intersecting and closer, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
          }
        }
      }
      walls = m_periodicWallTable[p1st->mySlot];
          // if it is NULL (was set in setupAfter... as initial value) we may skip it
      if(walls)
      {
        wallEnd = walls->end();
        for(wallIt = walls->begin(); wallIt != wallEnd; ++wallIt)
        {
//                 MSG_DEBUG("ValCalculatorSymmetryBCScalar::computeDists", "WALLISFIRST, periodic: checking, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
          if(wallIt->first->intersects
             (
                /*here we may always take [rj] since we used the periodic image of the 
             wall particle for the computation of its distance to the wall*/ 
             p2nd->r, 
                  /*const point_t&*/ pD->cartesian(), // !!! [rij] points from j to i !!!
                  /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
             )
            )
          {
            if(tempDist < innerDist)
            {
              innerDist = tempDist; 
              outerDist = wallIt->second;
	      m_currentWall = wallIt->first;
//                   MSG_DEBUG("ValCalculatorSymmetryBCScalar::computeDists", "WALLISFIRST, periodic: intersecting and closer, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
            }
          }
        }
      }
/*      else 
             MSG_DEBUG("ValCalculatorSymmetryBCScalar::computeDists", "WALLISFIRST, periodic: NO WALLS");*/
          // if all worked fine, we may now compute the value of the outer particle
          // the check for innerDist != HUGE_VAL is done below
          // the first term is for Dirichlet != 0. We assume the value was assigned 
          // to the wall particles.
//         pD->tag.pointByOffset(m_slot) =
//             (outerDist+innerDist)*(p1st->v)/innerDist-(outerDist/innerDist)*p2nd->v;
// //            MSG_DEBUG("ValCalculatorSymmetryBCScalar::compute", "WALLISFIRST: c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << ", tester=" << tester);

//     }
//     if(innerDist == HUGE_VAL)
//       throw gError("ValCalculatorSymmetryBCScalar::compute", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + "c1=" + ObjToString(pD->firstPart()->c) + ", c2=" + ObjToString(pD->secondPart()->c) + ", r1=" + ObjToString(pD->firstPart()->r) + ", r2=" + ObjToString(pD->secondPart()->r));
         
    //   } // this is the closing of the commented out if(m_cutoff > ...)                                                                      
    
    }
}


bool ValCalculatorSymmetryBCScalar::findStage()
{
  return Symbol::findStageNewPrelim();
}


bool ValCalculatorSymmetryBCScalar::findStage_0()
{
  return Symbol::findStageNewPrelim_0();
}
