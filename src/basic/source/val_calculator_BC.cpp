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


#include "val_calculator_BC.h"

#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "wall.h"
#include "pair_creator.h"
#include "cell.h"
// #include <utility>


#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculatorBC::ValCalculatorBC(string symbol)
  : ValCalculatorPair(symbol)
{

}

ValCalculatorBC::ValCalculatorBC(/*Node*/Simulation* parent)
  : ValCalculatorPair(parent)
{
  m_stage = 0;
  init();
}

void ValCalculatorBC::init()
{
  
  BOOLPC(frozen, m_frozen, "Are the BCs applied on free(\"0\", \"no\") or frozen(\"1\", \"yes\") boundary particles?");

  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the boundary value.");
  
  STRINGPC(wallSpecies, m_wallSpecies, 
           "Species of the wall particles."
          );

  m_symbolName = "BC";
  m_wallSpecies = "undefined";
  m_frozen = true;
  
#ifdef _OPENMP
  m_particleCalculator = false;
#endif
}

void ValCalculatorBC::setup()
{
  M_CONTROLLER->registerForSetupAfterParticleCreation(this);
  if(m_species.first == "undefined")
    throw gError("ValCalculatorBC::setup for " + m_properties.className(), "Attribute 'species1' has value \"undefined\" .");
  if(m_species.second == "undefined")
    throw gError("ValCalculatorBC::setup for " + m_properties.className(), "Attribute 'species2' has value \"undefined\" .");
  ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  // check whether wall particles are first or second in the pairs
  if(m_wallSpecies == cp->firstSpecies())
    m_wallIsSecond = false;
  else if(m_wallSpecies == cp->secondSpecies())
    m_wallIsSecond = true;
  else
    throw gError("ValCalculatorBC::setup for " + m_properties.className(), "Attribute wallSpecies = \"" + m_wallSpecies + "\" is neither equal to attribute species1 = \"" + m_species.first + "\" nor to species2 = \"" + m_species.second + "\".");

  m_wallColour = M_MANAGER->getColour(m_wallSpecies);

  //  m_datatype = DataFormat::POINT;
  
  ValCalculatorPair::setup();
  
}

void /*pair<size_t, size_t>*/ ValCalculatorBC::setSlot(ColourPair* cp, size_t& slot, bool oneProp)
{
  MSG_DEBUG("ValCalculatorBC::setSlot", "CALLED");
  m_slot = slot = cp->tagFormat().addAttribute
    ("ValCalculator_" + myName() + "_" + cp->toString(), m_datatype/*DataFormat::POINT*/, false, m_symbolName).offset;
}

// FIXME: inline ?
void ValCalculatorBC::computeDists(Pairdist* pD, double& innerDist, double& outerDist)
{
//   MSG_DEBUG("ValCalculatorBC::computeDists", "START");

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
    if(m_wallIsSecond)
    {
           // find the wall-segment where the rij-vector is hitting and the distance of ri to it
      wallEnd = m_wallTable[p2nd->mySlot]->end();
      for(wallIt = m_wallTable[p2nd->mySlot]->begin(); wallIt != wallEnd; ++wallIt)
      {
//  	MSG_DEBUG("ValCalculatorBC::computeDists", "checking a standard wall");
        if(wallIt->first->intersects
           (
	    /*this is either [ri] or if necessary its periodic image*/ 
	    p2nd->r + pD->cartesian(), 
            /*const point_t&*/ -pD->cartesian(), // !!! -[rij] points from i to j !!!
	    /*double&*/ tempDist
	    )
	   )
	  {
// 	    MSG_DEBUG("ValCalculatorBC::computeDists", "tempDist = " << tempDist);
	    if(tempDist < innerDist) {
	      innerDist = tempDist; 
	      outerDist = wallIt->second;
	    }
	  }
	// else MSG_DEBUG("ValCalculatorBC::computeDists", "no intersection with standard wall");
      }
      walls = m_periodicWallTable[p2nd->mySlot];
      // if it is NULL (was set in setupAfter... as initial value) we may skip it
//       MSG_DEBUG("ValCalculatorBC::computeDists", "periodic wall-array is " << walls);
      if(walls)
      {
//  	MSG_DEBUG("ValCalculatorBC::computeDists", "checking a periodic wall");
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
            }
          }
        }
      }

//       MSG_DEBUG("ValCalculatorBC::computeDists", "WALLISSECOND: c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << "r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist);

    }
    else // so !m_wallIsSecond
    {
           // find the wall-segment (-1)*rij-vector is hitting and the distance of rj to it
      wallEnd = m_wallTable[p1st->mySlot]->end();
      for(wallIt = m_wallTable[p1st->mySlot]->begin(); wallIt != wallEnd; ++wallIt)
      {
// 	MSG_DEBUG("ValCalculatorBC::computeDists", "WALLISFIRST, non-periodic: checking, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
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
// 	    MSG_DEBUG("ValCalculatorBC::computeDists", "WALLISFIRST, non-periodic: intersecting and closer, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
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
//                 MSG_DEBUG("ValCalculatorBC::computeDists", "WALLISFIRST, periodic: checking, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
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
//                   MSG_DEBUG("ValCalculatorBC::computeDists", "WALLISFIRST, periodic: intersecting and closer, c1=" << pD->firstPart()->c << ", c2=" << pD->secondPart()->c << ", r1=" << pD->firstPart()->r << ", r2=" << pD->secondPart()->r << "v1=" << pD->firstPart()->v << ", v2=" << pD->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
            }
          }
        }
      }
/*      else 
             MSG_DEBUG("ValCalculatorBC::computeDists", "WALLISFIRST, periodic: NO WALLS");*/

     if(innerDist == HUGE_VAL)
       throw gError("ValCalculatorBC::computeDists", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pD->firstPart()->mySlot) + ", slot2=" + ObjToString(pD->secondPart()->mySlot) + "c1=" + ObjToString(pD->firstPart()->c) + ", c2=" + ObjToString(pD->secondPart()->c) + ", r1=" + ObjToString(pD->firstPart()->r) + ", r2=" + ObjToString(pD->secondPart()->r));
         
    //   } // this is the closing of the commented out if(m_cutoff > ...)                                                                          
    }
}

// for each wall-particle, find the walls in range
void ValCalculatorBC::setupAfterParticleCreation()
{
  double range = M_PHASE->pairCreator()->interactionCutoff();
//   double range = M_SIMULATION->maxCutoff;

//   MSG_DEBUG("ValCalculatorDirichletBC::setupAfterParticleCreation", "START: m_wallColour = " << m_wallColour);
  assert(m_wallColour >= 0);
  Phase* phase = M_PHASE;
  ManagerCell* manager = M_MANAGER;
  cuboid_t box = phase->boundary()->boundingBox();
  point_t boxMiddle = 0.5*(box.corner1 + box.corner2);
  bool_point_t pFront = phase->boundary()->periodicityFront();
  bool_point_t pBack = phase->boundary()->periodicityBack();
  // reserve space for each particle
  // (faster than doing a push_back in the particle loop because allocation only once)
  if(m_frozen)
    {
      m_wallTable.resize(phase->returnNofFrozenPC(m_wallColour) + 
			 phase->returnNofFrozenPC(m_wallColour));
      // I have to reserve the memory here too because the slots must correspond to the particle slots
      m_periodicWallTable.resize(phase->returnNofFrozenPC(m_wallColour) + 
				 phase->returnNofFrozenPC(m_wallColour));
  if(!m_wallTable.size())
    throw gError("ValCalculatorBC::setupAfterParticleCreation", "No frozen particles of species " + m_wallSpecies + " found. Have you created them? Have you set attribute 'frozen' of this module (and possibly the ParticleCreator's) correctly?");    
    }
  else
    {
      m_wallTable.resize(phase->returnNofPartC(m_wallColour) + 
			 phase->returnNofPartC(m_wallColour));
      // I have to reserve the memory here too because the slots must correspond to the particle slots
      m_periodicWallTable.resize(phase->returnNofPartC(m_wallColour) + 
				 phase->returnNofPartC(m_wallColour));
  if(!m_wallTable.size())
    throw gError("ValCalculatorBC::setupAfterParticleCreation", "No free particles of species " + m_wallSpecies + " found. Have you created them? Have you set attribute 'frozen' of this module (and possibly the ParticleCreator's) correctly?");    
    }

  FOR_EACH_PARTICLE_C
      (phase,
       (size_t) m_wallColour,
       // allocate memory for the vector the current entry of m_wallTable (belonging to 
       // the current particle) is pointing to
       // we set the size to one since each particle should have at least one wall
       // allocation for more walls should happen so seldomly that it shouldn't matter
       m_wallTable[__iSLFE->mySlot] = NULL;
       // "()" necessary because the "," disturbs the compiler
       (m_wallTable[__iSLFE->mySlot] = new vector<pair<Wall*, double> >(0));
       assert(m_wallTable[__iSLFE->mySlot]->empty());
       // for the periodic table we create a vector<pair<...> > only when needed
       m_periodicWallTable[__iSLFE->mySlot] = NULL;
       point_t pos = __iSLFE->r;
       Cell* cell = manager->findCell(pos);
       if(!cell)
         throw gError("ValCalculatorBC::setupAfterParticleCreation", "FATAL ERROR: No cell found for boundary particle! Check your geometry and your settings. If this doesn't help, contact the programmers.\nParticle information:\nposition: " + ObjToString(__iSLFE->r) + "\nslot: " + ObjToString(__iSLFE->mySlot) );

       // loop over all walls of the cell and its neighbouring cells (not the periodic neighbours unfortunately!)
       list<Wall*>* walls = cell->allWalls();
       for(list<Wall*>::iterator wallIt = walls->begin(); wallIt != walls->end(); ++wallIt)
       {
         Wall* wall = *wallIt;
         if(wall->isInRange(range, pos))
         {
          // store a pointer to the wall and the distance of the particle to its plane
          // m_wallTable is vector<vector<pair>* >
           m_wallTable[__iSLFE->mySlot]->push_back( pair<Wall*, double>(wall, wall->distToPlane(pos)) );
         }
       } // loop over walls

       // for the periodic directions, we have to check if the periodic image is in range of the walls contained in the periodic neighbours of the current cell

       for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
	 if(pFront[dir]) {
	   if(!pBack[dir])
	     throw gError("ValCalculatorBC::setupAfterParticleCreation", "FATAL ERROR: Cannot handle different periodicities in the same direction! Check your geometry and your settings. If this doesn't help, contact the programmers.\nParticle information:\nposition: " + ObjToString(__iSLFE->r) + "\nslot: " + ObjToString(__iSLFE->mySlot) + "\ncurrent direction = " + ObjToString(dir)); 
	   // compute periodic image of pos, which is different for different 
	   // periodic directions
	   point_t posper = pos;
	   if(pos[dir] > boxMiddle[dir] ) {
	     // so particle is in the right half and must be shifted backwards 
	     posper[dir] -= (box.corner2[dir] - box.corner1[dir]);
	   }
	   else {
	     // so particle is in the left half and must be shifted forwards 
	     posper[dir] += (box.corner2[dir] - box.corner1[dir]);
	   }
	   
	   // additional "()"-brackets because compiler gets into trouble with the "," and the FOR_EACH_PARTICLE_C macro
	   // somehow I don't get this line compiled
	   // (list<pair<int, Cell* > > *indOutlets = cell->indirectOutlets());
	   // for a wall particle far away from periodic boundaries, 
	   // the cell shouldn't contain any indirect outlets, so nothing is done
	   for(list<pair< point_t, Cell*> >::iterator cellIt = cell->indirectOutlets()->begin(); cellIt != cell->indirectOutlets()->end(); ++cellIt) {

	     walls = cellIt->second->allWalls();
	     for(list<Wall*>::iterator wallIt = walls->begin(); wallIt != walls->end(); ++wallIt) {
	       Wall* wall = *wallIt;

	       if(wall->isInRange(range, posper)) {
		 // store a pointer to the wall and the distance of the particle to its plane		 
		 if(m_periodicWallTable[__iSLFE->mySlot]) { // already a saved wall?
		   // check if current wall is new
		   bool isNew = true;
		   for(vector<pair<Wall*, double> >::iterator wIt2 = 
			 m_periodicWallTable[__iSLFE->mySlot]->begin(); 
		       wIt2 != m_periodicWallTable[__iSLFE->mySlot]->end(); ++wIt2) {
		     if((*wall) == (*(wIt2->first))) {
		       isNew = false;
		     }
		   }
		   if(isNew) {
		     m_periodicWallTable[__iSLFE->mySlot]
		       ->push_back( pair<Wall*, double>(wall, wall->distToPlane(posper)) );
		   }
		 }
		 else {// so it is still empty => allocate and add wall 
		   (m_periodicWallTable[__iSLFE->mySlot] = new vector<pair<Wall*, double> >(0));
		 m_periodicWallTable[__iSLFE->mySlot]->push_back( pair<Wall*, double>(wall, wall->distToPlane(posper)) );		 
		 }
	       } // if(wall->isInRange(..))
	     } // loop over walls
	   } // loop over indirect outlets
	 } // if(pFront[dir]) 
       } // loop over SPACE_DIMS            

       if((*(m_wallTable[__iSLFE->mySlot])).empty())
         throw gError("ValCalculatorBC::setupAfterParticleCreation", "FATAL ERROR: Boundary particle without wall found! Check your geometry and your settings. If this doesn't help, contact the programmers.\nParticle information:\nposition: " + ObjToString(__iSLFE->r) + "\nslot: " + ObjToString(__iSLFE->mySlot) );

       );
}
