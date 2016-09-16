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

#include "linked_list_creator.h"
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

class ManagerCell;

/* Register this PairCreator with the factory. */
const PairCreator_Register<LinkedListCreator> linked_list_creator("LinkedListCreator");



#define M_PHASE  ((Phase*) m_parent)
#define M_MANAGER  M_PHASE->manager()


LinkedListCreator::LinkedListCreator(Phase* p): PairCreator(p)
{
  init();
}


LinkedListCreator::~LinkedListCreator()
{

}


void LinkedListCreator::init()
{
  m_properties.setClassName("LinkedListCreator");
  m_properties.setDescription("Creates linked list pairs.");

}

void LinkedListCreator::setup()
{
  PairCreator::setup();
}


void LinkedListCreator::invalidatePositions()
{
  // avoids pair-computation if there are no non-bonded pairs needed
  m_valid_dist = true;
  FOR_EACH_COLOUR_PAIR
  (M_MANAGER,
   // zeroing the pairs in the same loop
   for (size_t t = 0; t < global::n_threads; t++) {
     cp->freePairs()[t].clear();
     cp->frozenPairs()[t].clear();
   }
    m_valid_dist = !(cp->needPairs()) && m_valid_dist;
  );
}


void LinkedListCreator::createDistances()
{
    /* Can be called multiple times in one time step.
     Distances are updated only if invalidatePositions has been
     called. */
 if (!m_valid_dist) {

#ifdef _OPENMP
// int count =0;

#pragma omp for ordered
  for (int t = 0; t < global::n_threads; ++t) {
    CellLink* first = M_MANAGER->firstLink()[t];
    for (CellLink* cl = first; cl != NULL; cl = cl->next) {	   
      cl->createDistances(t);
//       ++count;
    }
  }
//  MSG_DEBUG("LinkedListCreator::createDistances", "count = " << count << " active links = " << M_MANAGER->activeLinks()[t]);

#else
  LL_FOR_EACH__PARALLEL
  (AbstractCellLink,
   M_MANAGER->firstLink(),
   M_MANAGER->activeLinks(),
   NULL,
   
     i->createDistances();
  );
     // MSG_DEBUG("LinkedListCreator::createDistances", "number of links = " << M_MANAGER->activeLinks());
#endif

    // randomize pairs if wished
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

    m_valid_dist = true;
  }
}
