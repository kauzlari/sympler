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


using namespace std;


#include "f_pair_vels_BC.h"

#include "threads.h"
#include "particle.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

const GenFTypeConcr<FPairVelsBC> f_pair_vels_BC("FPairVelsBC");

//---- Constructors/Destructor ----

FPairVelsBC::FPairVelsBC(Simulation *simulation): FPairArbitraryBC(simulation), m_wallColour(-1)
{
  init();
}


FPairVelsBC::~FPairVelsBC()
{
}


void FPairVelsBC::init()
{
  m_properties.setClassName("FPairVelsBC");

  m_properties.setDescription(
      "This is a completely general pair force FPV on a velocity v_i designed "
      "for boundary interactions including Dirichlet boundary conditions "
      "such that:\n"
      "dv_i = FPV*dt\n"
      "     = particleFactor_i(BV)*Sum_j(pairFactor_ij(BV)*weight_ij)*dt\n"
      "where pairFactor_ij(BV) includes all pair contributions of the pair ij,\n"
      "weight_ij represents -W'(rij)/rij of the used weighting function W.,\n"
      "particleFactor_i(j)(BV) represents factors specific to particle i(j),\n"
      "and BV is the symbol for the boundary value (usually of some "
      "boundary-particles), which may be used in the expressions as indicated.\n"
      "Notice that it is crucial to define the attribute 'wallSpecies' correctly. "
      "If no wall particles exist for the given species, the behaviour will be "
      "most likely unintentional. The same is true if non-wall particles of the "
      "same species exist."
  );

  m_symmetry = -1;

  STRINGPC(wallSpecies, m_wallSpecies,
          "Species of the wall particles."
         );
  m_wallSpecies = "undefined";
  m_pairFactorStr = "idVec(1)";
  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";
}


//---- Methods ----

#ifndef _OPENMP
void FPairVelsBC::computeForces(Pairdist* pair, int force_index)
#else
void FPairVelsBC::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{

//   M_PAIRCREATOR->createDistances();

  double innerDist;
  double outerDist;
  double tempDist;
  vector<std::pair<Wall*, double> >* walls;
  vector<std::pair<Wall*, double> >::iterator wallIt;
  vector<std::pair<Wall*, double> >::iterator wallEnd;
  
  if (this->m_cutoff > pair->abs())
    {
      point_t temp;
      
      point_t tester;
      double arg[3];
      
      Particle* p1st = pair->firstPart();
      Particle* p2nd = pair->secondPart();
      //          point_t hitpos;
      innerDist = HUGE_VAL;
      if(m_wallIsSecond)
	{
	  // find the wall-segment where the rij-vector is hitting and the distance of ri to it
	  //           MSG_DEBUG("FPairVelsBC::computeForces", "10");
          wallEnd = m_wallTable[p2nd->mySlot]->end();
          for(wallIt = m_wallTable[p2nd->mySlot]->begin(); wallIt != wallEnd; ++wallIt)
	    {
	      //             MSG_DEBUG("FPairVelsBC::computeForces", "20");
	      if(wallIt->first->intersects
		 (
		  /*this is either [ri] or if necessary its periodic image*/
		  p2nd->r + pair->cartesian(),
		  /*const point_t&*/ -pair->cartesian(), // !!! -[rij] points from i to j !!!
		  /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
		  )
		 )
		{
		  //               MSG_DEBUG("FPairVelsBC::computeForces", "30");
		  if(tempDist < innerDist)
		    {
		      innerDist = tempDist;
		      outerDist = wallIt->second;
		    }
		}
	    }
          walls = m_periodicWallTable[p2nd->mySlot];
          // if it is NULL (was set in setupAfter... as initial value) we may skip it
          if(walls)
	    {
	      wallEnd = walls->end();
	      //             MSG_DEBUG("FPairVelsBC::computeForces", "40");
	      for(wallIt = walls->begin(); wallIt != wallEnd; ++wallIt)
		{
		  //               MSG_DEBUG("FPairVelsBC::computeForces", "50");
		  if(wallIt->first->intersects
		     (
		      /*here we may always take [ri] since we used the periodic image of the
			wall particle for the computation of its distance to the wall*/
		      p1st->r,
		      /*const point_t&*/ -pair->cartesian(), // !!! -[rij] points from i to j !!!
		      /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
		      )
		     )
		    {
		      //                 MSG_DEBUG("FPairVelsBC::computeForces", "60");
		      if(tempDist < innerDist)
			{
			  innerDist = tempDist;
			  outerDist = wallIt->second;
			}
		    }
		}
	    }
          // if all worked fine, we may now compute the value of the outer particle
          // the check for innerDist != HUGE_VAL is done below
          for (size_t i = 0; i < 3; ++i)
	    {
	      //             MSG_DEBUG("FPairVelsBC::computeForces", "70");
	      // the first term is for Dirichlet != 0. We assume the value was assigned
	      // to the wall particles.
	      tester[i] = arg[i] = (outerDist+innerDist)*(p2nd->v[i])/innerDist-(outerDist/innerDist)*p1st->v[i];
	    }
	  
	  //           MSG_DEBUG("FPairVelsBC::computeForces", "WALLISSECOND: c1=" << pair->firstPart()->c << ", c2=" << pair->secondPart()->c << "r1=" << pair->firstPart()->r << ", r2=" << pair->secondPart()->r << "v1=" << pair->firstPart()->v << ", v2=" << pair->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << ", tester=" << tester);
	  
	}
      else // so !m_wallIsSecond
	{
	  // find the wall-segment (-1)*rij-vector is hitting and the distance of rj to it
	  wallEnd = m_wallTable[p1st->mySlot]->end();
	  for(wallIt = m_wallTable[p1st->mySlot]->begin(); wallIt != wallEnd; ++wallIt)
	    {
	      //              MSG_DEBUG("FPairVelsBC::computeForces", "WALLISFIRST, non-periodic: checking, c1=" << pair->firstPart()->c << ", c2=" << pair->secondPart()->c << ", r1=" << pair->firstPart()->r << ", r2=" << pair->secondPart()->r << "v1=" << pair->firstPart()->v << ", v2=" << pair->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
	      if(wallIt->first->intersects
		 (
		  /*this is either [rj] or if necessary its periodic image*/
		  p1st->r - pair->cartesian(),
		  /*const point_t&*/ pair->cartesian(), // !!! [rij] points from j to i !!!
		  /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
		  )
		 )
		{
		  if(tempDist < innerDist)
		    {
		      innerDist = tempDist;
		      outerDist = wallIt->second;
		      //                  MSG_DEBUG("FPairVelsBC::computeForces", "WALLISFIRST, non-periodic: intersecting and closer, c1=" << pair->firstPart()->c << ", c2=" << pair->secondPart()->c << ", r1=" << pair->firstPart()->r << ", r2=" << pair->secondPart()->r << "v1=" << pair->firstPart()->v << ", v2=" << pair->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
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
		  //                MSG_DEBUG("FPairVelsBC::computeForces", "WALLISFIRST, periodic: checking, c1=" << pair->firstPart()->c << ", c2=" << pair->secondPart()->c << ", r1=" << pair->firstPart()->r << ", r2=" << pair->secondPart()->r << "v1=" << pair->firstPart()->v << ", v2=" << pair->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
		  if(wallIt->first->intersects
		     (
		      /*here we may always take [rj] since we used the periodic image of the
			wall particle for the computation of its distance to the wall*/
		      p2nd->r,
		      /*const point_t&*/ pair->cartesian(), // !!! [rij] points from j to i !!!
		      /*point_t&*/ /*hit_pos,*/ /*double&*/ tempDist
		      )
		     )
		    {
		      if(tempDist < innerDist)
			{
			  innerDist = tempDist;
			  outerDist = wallIt->second;
			  //                   MSG_DEBUG("FPairVelsBC::computeForces", "WALLISFIRST, periodic: intersecting and closer, c1=" << pair->firstPart()->c << ", c2=" << pair->secondPart()->c << ", r1=" << pair->firstPart()->r << ", r2=" << pair->secondPart()->r << "v1=" << pair->firstPart()->v << ", v2=" << pair->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << "\nwall: " << wallIt->first->toString());
			}
		    }
		}
	    }
	  //            else
	  //              MSG_DEBUG("FPairVelsBC::computeForces", "WALLISFIRST, periodic: NO WALLS");
          // if all worked fine, we may now compute the value of the outer particle
          // the check for innerDist != HUGE_VAL is done below
	  for (size_t i = 0; i < 3; ++i)
	    {
	      // the first term is for Dirichlet != 0. We assume the value was assigned
	      // to the wall particles.
	      tester[i] = arg[i] = (outerDist+innerDist)*(p1st->v[i])/innerDist-(outerDist/innerDist)*p2nd->v[i];
	    }
	  //            MSG_DEBUG("FPairVelsBC::computeForces", "WALLISFIRST: c1=" << pair->firstPart()->c << ", c2=" << pair->secondPart()->c << ", r1=" << pair->firstPart()->r << ", r2=" << pair->secondPart()->r << "v1=" << pair->firstPart()->v << ", v2=" << pair->secondPart()->v << ", iDist=" << innerDist << ", oDist=" << outerDist << ", tester=" << tester);
	  
	}
      if(innerDist == HUGE_VAL)
	throw gError("FPairVelsBC::computeForces", "No wall found for pair. Check your geometry and other settings. If this doesn't help, contact the programmers. \nDetails: slot1=" + ObjToString(pair->firstPart()->mySlot) + ", slot2=" + ObjToString(pair->secondPart()->mySlot) + "c1=" + ObjToString(pair->firstPart()->c) + ", c2=" + ObjToString(pair->secondPart()->c) + ", r1=" + ObjToString(pair->firstPart()->r) + ", r2=" + ObjToString(pair->secondPart()->r));
            
    this->m_pairFactor(&temp, arg, &(*pair));

    point_t fi;
    point_t fj;

    // compute the particle-expressions
    this->m_1stparticleFactor(&fi, arg, &(*pair));
    this->m_2ndparticleFactor(&fj, arg, &(*pair));

    // loop necessary because operator* of math_vector_t does scalar product
    for(size_t i = 0; i < SPACE_DIMS; ++i)
    {
      fi[i] *= temp[i];
      fj[i] *= temp[i];
    }

    fi *= this->m_wf->weight(pair, pair->secondPart()->r);
    fj *= this->m_wf->weight(pair, pair->firstPart()->r);

#ifndef _OPENMP
    if (pair->actsOnFirst())
      pair->firstPart()->force[/*self->m_*/force_index] += fi;
    if (pair->actsOnSecond())
      pair->secondPart()->force[/*self->m_*/force_index] += m_symmetry*fj;
#else
    if (pair->actsOnFirst()) {
      for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
        (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi[_i];
      }
    }
    if (pair->actsOnSecond()) {
      for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
        (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += m_symmetry*fj[_i];
      }
    }
#endif

       }
}


void FPairVelsBC::computeForces(Particle* part, int force_index)
{
  throw gError("FPairVelsBC::computeForces", "Fatal error: do not call FPairVelsBC::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairVelsBC::computeForces(int force_index)
{
  throw gError("FPairVelsBC::computeForces", "Fatal error: do not call FPairVelsBC::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairVelsBC::setup()
{
  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  M_CONTROLLER->registerForSetupAfterParticleCreation(this);
  FPairArbitraryBC::setup();

  m_pairFactor.setReturnType(Variant::VECTOR);
  m_pairFactor.addVariable("[BV]");
  m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_1stparticleFactor.addVariable("[BV]");
  m_2ndparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.addVariable("[BV]");

  // check whether wall particles are first or second in the pairs
  if(m_wallSpecies == m_cp->firstSpecies())
    m_wallIsSecond = false;
  else if(m_wallSpecies == m_cp->secondSpecies())
    m_wallIsSecond = true;
  else
    throw gError("FPairVelsBC::setup", "Attribute wallSpecies = \"" + m_wallSpecies + "\" is neither equal to attribute species1 = \"" + m_species.first + "\" nor to species2 = \"" + m_species.second + "\".");

  m_wallColour = M_MANAGER->getColour(m_wallSpecies);


// Setting the colours this force works on
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if ((col1 == cp->secondColour()) && (col2 == cp->firstColour())) {
    size_t dummy = col1;
    col1 = col2;
    col2 = dummy;
  }
  else if ((col1 == cp->firstColour()) && (col2 == cp->secondColour())) {
//    doesn't change
  }
  else {
    throw gError("FPairVelsBC::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FPairVelsBC::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof = "vel_pos";
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
      if (col1 == cp->firstColour()) {
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();
      }
      else if (col1 == cp->secondColour()) {
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();      
      }
      else {
        throw gError("FPairVelsBC::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
  if (col2 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
      if (col2 == cp->secondColour()) {
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();
      }
      else if (col2 == cp->firstColour()) {
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();      
      }
      else {
        throw gError("FPairVelsBC::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
}
#endif


// for each wall-particle, find the walls in range
void FPairVelsBC::setupAfterParticleCreation()
{
  MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "START: m_wallColour = " << m_wallColour);
  assert(m_wallColour >= 0);
  Phase* phase = M_PHASE;
  ManagerCell* manager = M_MANAGER;
  cuboid_t box = phase->boundary()->boundingBox();
  point_t boxMiddle = 0.5*(box.corner1 + box.corner2);
  bool_point_t pFront = phase->boundary()->periodicityFront();
  bool_point_t pBack = phase->boundary()->periodicityBack();
  // reserve space for each particle
  // (faster than doing a push_back in the particle loop because allocation only once)
  m_wallTable.resize(phase->returnNofFrozenPC(m_wallColour) +
      phase->returnNofFrozenPC(m_wallColour));
  // I have to reserve the memory here too because the slots must correspond to the particle slots
  m_periodicWallTable.resize(phase->returnNofFrozenPC(m_wallColour) +
      phase->returnNofFrozenPC(m_wallColour));
  if(!m_wallTable.size())
    throw gError("FPairVelsBC::setupAfterParticleCreation", "No particles of species " + m_wallSpecies + " found.");
  FOR_EACH_PARTICLE_C
      (phase,
       (size_t) m_wallColour,
       // allocate memory for the vector the current entry of m_wallTable (belonging to
       // the current particle) is pointing to
       // we set the size to one since each particle should have at least one wall
       // allocation for more walls should happen so seldomly that it shouldn't matter
       m_wallTable[__iSLFE->mySlot] = NULL;
       (m_wallTable[__iSLFE->mySlot] = new vector<pair<Wall*, double> >(0));
       assert(m_wallTable[__iSLFE->mySlot]->empty());
       // for the periodic table we create a vector<pair<...> > only when needed
       m_periodicWallTable[__iSLFE->mySlot] = NULL;
       point_t pos = __iSLFE->r;
       Cell* cell = manager->findCell(pos);
       if(!cell)
         throw gError("FPairVelsBC::setupAfterParticleCreation", "FATAL ERROR: No cell found for boundary particle! Check your geometry and your settings. If this doesn't help, contact the programmers.\nParticle information:\nposition: " + ObjToString(__iSLFE->r) + "\nslot: " + ObjToString(__iSLFE->mySlot) );
       // loop over all walls of the cell and its neighbouring cells
       list<Wall*>* walls = cell->allWalls();
       for(list<Wall*>::iterator wallIt = walls->begin(); wallIt != walls->end(); ++wallIt)
       {
         Wall* wall = *wallIt;
//          MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Before ifinrange");
         if(wall->isInRange(m_cutoff, pos))
         {
          // store a pointer to the wall and the distance of the particle to its plane
          // m_wallTable is vector<vector<pair>* >
           m_wallTable[__iSLFE->mySlot]->push_back( pair<Wall*, double>(wall, wall->distToPlane(pos)) );
//            MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Found wall: p=" + ObjToString(__iSLFE->r) + ", pslot=" + ObjToString(__iSLFE->mySlot) + ", wall=" + wall->toString() + ", wallTable-size = " + ObjToString(m_wallTable[__iSLFE->mySlot]->size()));
         }
         else
         {
           // for the periodic directions, we have to check if the periodic image is in range
//            MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Before complicated periodic debug");
//            MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Searching periodic wall: p=" + ObjToString(__iSLFE->r) + ", pslot=" + ObjToString(__iSLFE->mySlot) + ", wall=" + wall->toString() + ", periodicWallTable-pointer = " + ObjToString(m_periodicWallTable[__iSLFE->mySlot]));

/*           if(m_periodicWallTable[__iSLFE->mySlot])
             MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Searching periodic wall, periodicWallTable-size = " +
            ObjToString(m_periodicWallTable[__iSLFE->mySlot]->size()));
*/
           for(size_t dir = 0; dir < SPACE_DIMS; ++dir)
           {
             if(pFront[dir])
             {
//                MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Searching periodic wall: dir=" + ObjToString(dir) + " is periodic.");
               if(!pBack[dir])
                 throw gError("FPairVelsBC::setupAfterParticleCreation", "FATAL ERROR: Cannot handle different periodicities in the same direction! Check your geometry and your settings. If this doesn't help, contact the programmers.\nParticle information:\nposition: " + ObjToString(__iSLFE->r) + "\nslot: " + ObjToString(__iSLFE->mySlot) + "\nwall=" + wall->toString() + "\ncurrent direction = " + ObjToString(dir));
              // compute periodic image of pos, which is different for different
              // periodic directions
              point_t posper = pos;
              if(pos[dir] > boxMiddle[dir] )
              {
                // so particle is in the right half and must be shifted backwards
                posper[dir] -= (box.corner2[dir] - box.corner1[dir]);
              }
              else
              {
                // so particle is in the left half and must be shifted forwards
                posper[dir] += (box.corner2[dir] - box.corner1[dir]);
              }
//               MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Searching periodic wall: posper is now " + ObjToString(posper));
              // is the periodic image in range of the wall?
              if(wall->isInRange(m_cutoff, posper))
              {
//                 MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Searching periodic wall: posper is in range");
                if(m_periodicWallTable[__iSLFE->mySlot]) // already a saved wall?
                {
                  // check if current wall is new
                  bool isNew = true;
                  for(vector<pair<Wall*, double> >::iterator wIt2 =
                      m_periodicWallTable[__iSLFE->mySlot]->begin();
                      wIt2 != m_periodicWallTable[__iSLFE->mySlot]->end(); ++wIt2)
                  {
                    if((*wall) == (*(wIt2->first)))
                    {
                      isNew = false;
//                       MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Found already existing wall; NOT ADDING: p=" + ObjToString(__iSLFE->r) + ", pslot=" + ObjToString(__iSLFE->mySlot) + ", wall=" + wall->toString() + ", periodicWallTable-size = " + ObjToString(m_periodicWallTable[__iSLFE->mySlot]->size()));
                    }
                  }
                  if(isNew)
                  {
                    m_periodicWallTable[__iSLFE->mySlot]
                        ->push_back( pair<Wall*, double>(wall, wall->distToPlane(posper)) );
//                     MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Found OTHER periodic wall: p=" + ObjToString(__iSLFE->r) + ", pslot=" + ObjToString(__iSLFE->mySlot) + ", wall=" + wall->toString() + ", periodicWallTable-size = " + ObjToString(m_periodicWallTable[__iSLFE->mySlot]->size()));
                  }
                }
                else // so it is still empty => allocate and create
                {
                  (m_periodicWallTable[__iSLFE->mySlot] = new vector<pair<Wall*, double> >(0));
                  m_periodicWallTable[__iSLFE->mySlot]
                      ->push_back( pair<Wall*, double>(wall, wall->distToPlane(posper)) );
//                   MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Found FIRST periodic wall: p=" + ObjToString(__iSLFE->r) + ", pslot=" + ObjToString(__iSLFE->mySlot) + ", wall=" + wall->toString() + ", periodicWallTable-size = " + ObjToString(m_periodicWallTable[__iSLFE->mySlot]->size()));

                }
              }
//               else
//                 MSG_DEBUG("FPairVelsBC::setupAfterParticleCreation", "Searching periodic wall: posper is NOT in range");

             }
           }
         }
       }
       if((*(m_wallTable[__iSLFE->mySlot])).empty())
          throw gError("FPairVelsBC::setupAfterParticleCreation", "FATAL ERROR: Boundary particle without wall found! Check your geometry and your settings. If this doesn't help, contact the programmers.\nParticle information:\nposition: " + ObjToString(__iSLFE->r) + "\nslot: " + ObjToString(__iSLFE->mySlot) );
      );
}
