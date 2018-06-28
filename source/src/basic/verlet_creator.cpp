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
#include <limits>

#include "verlet_creator.h"
#include "threads.h"
#include "inlet_cell.h"
#include "simulation.h"
#include "vertex_list.h"
#include "manager_cell.h"
#include "particle_cache.h"
#include "particle_creator.h"
#include "colour_pair.h"

#ifdef _OPENMP
  #include "omp.h"
#endif

#include "time.h"
// necessary for arguments that are pointers to Phase
#include "phase.h"


using namespace std;

#define M_PHASE  ((Phase*) m_parent)
#define M_MANAGER  M_PHASE->manager()
#define M_SIMULATION ((Simulation*) M_PHASE->parent())
#define M_CONTROLLER  M_SIMULATION->controller()

class ManagerCell;

/* Register this PairCreator with the factory. */
const PairCreator_Register<VerletCreator> verlet_creator("VerletCreator");


VerletCreator::VerletCreator(Phase* p): PairCreator(p)
{
  init();
}


VerletCreator::~VerletCreator()
{

}


void VerletCreator::init()
{
  m_properties.setClassName("VerletCreator");
  m_properties.setDescription
    ("Implementation of a PairCreator creating a neighbour list by "
     "applying first a cell subdivision and then the standard Verlet "
     "neighbour list method.\n"
     "Please NOTE: This PairCreator does not reset (to zero) symbols "
     "stored in "
     "pairs unless the whole neighbour list needs to be rebuilt. "
     "Therefore, while the distance vector [rij] of the pairs may "
     "already have been updated, the pairs may still hold data for a "
     "while (until the next call of Symbols responsible for the "
     "data), which was computed at an older pair state corresponding "
     "to an older [rij]. If the data depends on [rij] this may be "
     "undesirable, otherwise probably not. In any case, the user "
     "should consult the output of 'sympler --help workflow' to "
     "figure out, which symbols should ideally be computed at which "
     "stage. Note that data stored in the pairs may therefore NOT be "
     "identical at all instances during a time step to those from a "
     "simulation using a "
     "different PairCreator. Unthoughtful usage of the Symbol stages "
     "may therefore lead to PairCreator-dependent (and usually wrong) "
     "simulation results."
);

  DOUBLEPC
    (skinSize,
     m_skin_size,
     -HUGE_VAL,
     "The size of the additional skin added to the cutoff.");

  INTPC
    (every, m_every, -1, "If this value is larger than zero, the "
     "neighbour list is updated strictly every 'every' time steps. "
     "This update strategy requires some previous knowledge, "
     "otherwise you risk to miss neighbours if 'every' is too large.");

  STRINGPC
    (displacement, m_displacement_name,
     "Name of the particle displacement required to determine when "
     "the neighbour list must be updated. This quantity must be "
     "computed externally for each particle, for example by a "
     "suitable Integrator such as IntegratorVelocityVerletDisp.");

  m_displacement_name = "undefined";
  m_skin_size = -1;
  m_every = 0;

  // needed if m_every > 0
  m_counter = 0;

}

void VerletCreator::setup()
{
  PairCreator::setup();
  
  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);
  
  if (m_skin_size < 0.)
    throw gError("VerletCreator::setup", "Please define a positive skin size with attribute 'skinSize'!");

  if (m_every < 0)
    throw gError("VerletCreator::setup", "Please vive a non-negative value for attribute 'every'! Value \"" + ObjToString(m_every) + "\" is not allowed!");
  
  // we must set the cutoff before the particle creation and we can do it here because this setup is called after the setup of Forces, Symbols, and(!) Callables (e.g., Thermostats).
  FOR_EACH_COLOUR_PAIR
    (M_MANAGER,
     if(cp->needPairs())     
       cp->setCutoff(cp->cutoff() + m_skin_size);
     MSG_DEBUG("VerletCreator::setup", "CP(" << cp->firstSpecies() << "," << cp->secondSpecies() << "): now cp->cutoff()=" << cp->cutoff() << ", M_SIMULATION->maxCutoff=" << M_SIMULATION->maxCutoff); 
     );

 
  if(m_displacement_name == "undefined")
    throw gError("VerletCreator::setup", "Attribute 'displacement' has value \"undefined\". Please give the name under which the displacement will be computed by another module, e.g., by an IntegratorVelocityVerletDisp.");

  // set array sizes for simplicity to the number of colours. But we 
  // will not create the variables for all colours !!!  
  m_needDisp.resize(M_MANAGER->nColours());
  m_displacement_o.resize(M_MANAGER->nColours());
  m_displacementOld_o.resize(M_MANAGER->nColours());

  string displacementOld_name = string(m_displacement_name + "__Old");
  
  for (size_t colour = 0; colour < M_MANAGER->nColours(); ++colour) {
    
    // check which colours have an IntegratorPosition and hence need a displacement to be checked
    if(M_CONTROLLER->findIntegrator("IntegratorPosition", M_MANAGER->species(colour))) {
      MSG_DEBUG("VerletCreator::setup", "IntegratorPosition found for colour " << colour);
      m_needDisp[colour] = true;
      if (Particle::s_tag_format[colour].attrExists(m_displacement_name)) {
	DataFormat::attribute_t attrDisp = Particle::s_tag_format[colour].attrByName(m_displacement_name);
	
	if(attrDisp.datatype != DataFormat::POINT)
	  throw gError("VerletCreator::setup", "the symbol " + m_displacement_name + " is not registered as a point");
	
	m_displacement_o[colour] = Particle::s_tag_format[colour].attrByName(m_displacement_name).offset;
      }
      else {
	throw gError("VerletCreator::setup", "No symbol '" + m_displacement_name + "' found for species '" + M_MANAGER->species(colour) + "'! You need another module introducing it. This might be for example an IntegratorVelocityVerletDisp.");
      }
      
      if (!Particle::s_tag_format[colour].attrExists(displacementOld_name)) {      
	m_displacementOld_o[colour] = 
	  Particle::s_tag_format[colour].addAttribute
	  (displacementOld_name,
	   DataFormat::POINT,
	   true, // must be persistent!
	   displacementOld_name).offset;
      }
      else {
	throw gError("VerletCreator::setup", "Symbol '" + displacementOld_name + "' already existing for species '" + M_MANAGER->species(colour) + "' but I have to register mine! You may not define it twice, so please change the other symbol's name!");
      }
      
    }
    else {// no IntegratorPosition
      // dummy value that should never be used
      m_displacement_o[colour] = std::numeric_limits<size_t>::max();
      m_displacementOld_o[colour] = std::numeric_limits<size_t>::max();
      m_needDisp[colour] = false;
    }
       
  } // end for (size_t colour = 0; colour < M_MANAGER->nColours(); ++colour)
}

void VerletCreator::setupAfterParticleCreation()
{

  // initialise old displacement to 2*m_skinSize so that a neighbour-list 
  // is constructed in the first time step
  for (size_t col = 0; col < M_MANAGER->nColours(); ++col) {
    // check important because for those colours where not needed, the variable was not created 
    if(m_needDisp[col]) {
      FOR_EACH_FREE_PARTICLE_C__PARALLEL
	(M_PHASE, col, this,
	 point_t& oldDisp = i->tag.pointByOffset(m_displacementOld_o[col]);
	 for(size_t dir = 0; dir < SPACE_DIMS; ++dir)
	   oldDisp[dir] = 2*m_skin_size; 
	 );
    }
  }
}


void VerletCreator::invalidatePositions()
{
  // avoids pair-computation if there are no non-bonded pairs needed
  m_valid_dist = true;
  FOR_EACH_COLOUR_PAIR
  (M_MANAGER,
    m_valid_dist = !(cp->needPairs()) && m_valid_dist;
  );
}

void VerletCreator::createDistances()
{
  /* Can be called multiple times in one time step.
     Distances are recomputed always (unless m_valid_dist is true), 
     but the pair-list is rebuilt only if particle displacments 
     are too large */
  if (!m_valid_dist) {

    bool newList = false;

    if(m_every) {
      // we strictly update every m_every
      newList = !m_counter || (m_counter == size_t(m_every));
      // the following commented-out code is the long form of the one-liner
//       if(!m_counter)
// 	// so we are at the beginning and need a fresh pair list
// 	newList = true;
//       else {
// 	newList = (m_counter == m_every);
//       }
//       MSG_DEBUG("VerletCreator::createDistances", "m_every-case: newList=" << newList << "m_counter=" << m_counter << ", m_every=" << m_every);
    }
    else {
      // we update if the two largest displacements say so
      size_t nColours = M_MANAGER->nColours();
      Phase *phase = M_PHASE;
      
      double max_disp = 0;
      double max2 = 0;
      
      // The initial "old" displacements are equal to 2*m_skinSize
      // So, initially the neighbour-list will be constructed      
      for (size_t c = 0; c < nColours; ++c) {
	if(m_needDisp[c]) {
	  FOR_EACH_FREE_PARTICLE_C__PARALLEL
	    (phase, c, this,
	     point_t dispNow = (i->tag.pointByOffset(this->m_displacement_o[c]))-(i->tag.pointByOffset(this->m_displacementOld_o[c]));
	     double tempDisp = dispNow.abs();
	     
	     // 	 MSG_DEBUG("VerletCreator::createDistances", "test-output for testing break-command: at particle: " << __iSLFE->mySlot << ", max_disp_old=" << max_disp << ", max2_old=" << max2 << ", tempDisp="<< tempDisp);
	     
	     if (max_disp < tempDisp) {
	       // 	   MSG_DEBUG("VerletCreator::createDistances", "max_disp < tempDisp, max_disp_old=" << max_disp << ", tempDisp="<< tempDisp);
	       max_disp = tempDisp;
	     }
	     else if (max2 < tempDisp) {
	       // 	   MSG_DEBUG("VerletCreator::createDistances", "max2 < tempDisp, max2_old=" << max2 << ", tempDisp="<< tempDisp);
	       max2 = tempDisp;       
	     }
	     newList = (max_disp + max2) >= m_skin_size;
	     // breaks particle loop
	     if(newList) break;      
	     );
	} // end of if(m_needDisp[c]
	// for breaking colour loop
	if(newList) {
	  break;
	}      
      } // end of for (size_t c = 0; c < nColours; ++c)
    }

    if(newList) {// new neighbour-lists have to be created

//       MSG_DEBUG("VerletCreator::createDistances", "creating new neighbour list after " << m_counter << " steps.");
      m_counter = 1;

      ManagerCell* manager = M_MANAGER;

      size_t nColours = manager->nColours();
      Phase *phase = M_PHASE;

      // set old displacmements to new values
      for (size_t c = 0; c < nColours; ++c) {
	// check important because for those colours where not needed, the variable does not exist
	if(m_needDisp[c]) {
	  FOR_EACH_FREE_PARTICLE_C__PARALLEL
	    (phase, c, this,
	     (i->tag.pointByOffset(this->m_displacementOld_o[c])) = (i->tag.pointByOffset(this->m_displacement_o[c]));
	     );
	}
      }

      // clear old neighbour-lists
      FOR_EACH_COLOUR_PAIR
	(manager,
	 for (int t = 0; t < global::n_threads; ++t) {
	   cp->freePairs()[t].clear();
	   cp->frozenPairs()[t].clear();
	 }
	 );

      // create new neighbour list    
#ifdef _OPENMP
      // int count =0;
      
#pragma omp for ordered
      // FIXME!: so don't we distinguish between active and non-active links in the parallel case?
      for (int t = 0; t < global::n_threads; ++t) {
	CellLink* first = manager->firstLink()[t];
	for (CellLink* cl = first; cl != NULL; cl = cl->next) {	   
	  cl->createDistances(t);
	  //       ++count;
	}
      }
      //  MSG_DEBUG("VerletCreator::createDistances", "count = " << count << " active links = " << M_MANAGER->activeLinks()[t]);
      
#else
      LL_FOR_EACH__PARALLEL
	(CellLink,
	 manager->firstLink(),
	 manager->activeLinks(),
	 NULL,
	 
	 i->createDistances();
	 );
      // MSG_DEBUG("VerletCreator::createDistances", "number of links = " << M_MANAGER->activeLinks());
#endif
      
    } // end of if(newList)

    else {
      // just compute new distances of existing pairs

      ++m_counter;

      ManagerCell* manager = M_MANAGER;

      point_t boxSize = M_PHASE->boundary()->boundingBox().size();
      double size;

#pragma omp parallel for
      for (int t = 0; t < global::n_threads; ++t) {
	vector<ColourPair*>::iterator __end = manager->colourPairs().end();   
	for(vector<ColourPair*>::iterator __cp = manager->colourPairs().begin(); __cp != __end; ++__cp) {                                                
	  ColourPair *cp = *__cp; 
	  
	  if(cp->freePairsRandom(t).size())	{
	    for (PrimitiveSLEntry<size_t> *i = cp->freePairsRandom(t).first(); i != NULL; i = i->next) {
	      Pairdist* pair = &(cp->freePairs()[t][i->m_val]);
	      
	      pair->calculateCartDistance();
	      // periodic BCs
	      point_t& cartesian = pair->m_distance.cartesian;
	      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
		size = boxSize[dir];

		if(cartesian[dir] > 0.5*size) cartesian[dir] -= size; 
		if(cartesian[dir] < -0.5*size) cartesian[dir] += size; 
	      }
	      // absolute values
	      pair->m_distance.calcAbs();

	    }
	  }
	  else {
	    for (Pairdist *i = cp->freePairs()[t].first(); i != NULL; i = i->next) {
	      
	      Pairdist* pair = i;
	      
	      pair->calculateCartDistance();
	      // periodic BCs
	      point_t& cartesian = pair->m_distance.cartesian;
	      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
		size = boxSize[dir];

		if(cartesian[dir] > 0.5*size) cartesian[dir] -= size; 
		if(cartesian[dir] < -0.5*size) cartesian[dir] += size; 
	      }
	      // absolute values
	      pair->m_distance.calcAbs();

	    } 	  
	  }
	  
	  if(cp->frozenPairsRandom(t).size()) {
	    for (PrimitiveSLEntry<size_t> *i = cp->frozenPairsRandom(t).first(); i != NULL; i = i->next) {
	      Pairdist* pair = &(cp->frozenPairs()[t][i->m_val]);
	      
	      pair->calculateCartDistance();
	      // periodic BCs
	      point_t& cartesian = pair->m_distance.cartesian;
	      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
		size = boxSize[dir];

		if(cartesian[dir] > 0.5*size) cartesian[dir] -= size; 
		if(cartesian[dir] < -0.5*size) cartesian[dir] += size; 
	      }
	      // absolute values
	      pair->m_distance.calcAbs();

	    }
	  }
	  else {
	    for (Pairdist *i = cp->frozenPairs()[t].first(); i != NULL; i = i->next) {
	      Pairdist* pair = i;
	      
	      pair->calculateCartDistance();
	      // periodic BCs
	      point_t& cartesian = pair->m_distance.cartesian;
	      for(size_t dir = 0; dir < SPACE_DIMS; ++dir) {
		size = boxSize[dir];

		if(cartesian[dir] > 0.5*size) cartesian[dir] -= size; 
		if(cartesian[dir] < -0.5*size) cartesian[dir] += size; 
	      }
	      // absolute values
	      pair->m_distance.calcAbs();

	    }
	  }  
	}
      } // end of for (int t = 0; t < global::n_threads; ++t)
      
      
    } // end of else of if(newList)
    
    m_valid_dist = true;
  } // end of if(!m_valid_dist)
  
  
  // randomize pairs if wished (also if no new neighbour list has been created)
  if(M_PHASE->randomPairs())
    {
      FOR_EACH_COLOUR_PAIR
	(M_MANAGER,
	 
	 for (int t = 0; t < global::n_threads; ++t) {
	   
	   // free pairs
	   size_t listSize = cp->freePairs()[t].size();
	   if(listSize)
	     {
	       cp->freePairsRandom(t).clear();
	       for(size_t __i = 0; __i < listSize; ++__i)
		 cp->freePairsRandom(t).newEntry().m_val = __i;
	       for(size_t __i = listSize-1; __i > 0; --__i)
		 {
		   size_t slot = size_t(M_PHASE->rng().uniform()*(__i+1));
		   size_t temp = cp->freePairsRandom(t)[__i].m_val;
		   cp->freePairsRandom(t)[__i].m_val = cp->freePairsRandom(t)[slot].m_val;
		   cp->freePairsRandom(t)[slot].m_val = temp;
		 }
	     }
	   // frozen pairs
	   listSize = cp->frozenPairs()[t].size();
	   if(listSize)
	     {
	       cp->frozenPairsRandom(t).clear();
	       for(size_t __i = 0; __i < listSize; ++__i)
		 cp->frozenPairsRandom(t).newEntry().m_val = __i;
	       for(size_t __i = listSize-1; __i > 0; --__i)
		 {
		   size_t slot = size_t(M_PHASE->rng().uniform()*(__i+1));
		   size_t temp = cp->frozenPairsRandom(t)[__i].m_val;
		   cp->frozenPairsRandom(t)[__i].m_val = cp->frozenPairsRandom(t)[slot].m_val;
		   cp->frozenPairsRandom(t)[slot].m_val = temp;
		 }
	     }
	 }
	 );
    }
  
}
