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



#ifndef __VERLET_CREATOR_H
#define __VERLET_CREATOR_H

using namespace std;


#include <list>

#include "cell.h"
#include "gen_f.h"
#include "pairdist.h"
#include "pair_list.h"
#include "wall_container.h"
#include "pair_creator.h"
#include "node_many_children.h"


/*!
 * Implementation of a \a PairCreator creating a \a PairList by applying 
 * first a cell subdivision and then the standard Verlet-list method 
 */
class VerletCreator: public PairCreator
{
 protected:
  
  /*!
   * The skin size that is added to the cut-off
   */
  double m_skin_size;
  
  /*!
   * Number of time-steps after which the neighbour list should be updated. 
   * If zero, then the two largest particle-displacements decide whether 
   * there should be an update.
   */
  int m_every;
  
  /*!
   * The name of the displacement
   */
  string m_displacement_name;
  
  /*!
   * Do we need to check displacements for this colour?
   */
  vector<bool> m_needDisp;

  /*!
   * Tag offset of the displacement per colour
   */
  vector<size_t> m_displacement_o;

  /*!
   * tag offset of helper displacement per colour stored at that time the
   * \a neighbour list was identified to require an update
   */
  vector<size_t> m_displacementOld_o;

  /*!
   * Initialize
   */
  void init();

 private:

  /*!
   * Helper for tracking and debugging purposes
   */
  size_t m_counter;

 public:

  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   */
  VerletCreator(Phase* p);

  /*!
   * Destructor
   */
  virtual ~VerletCreator();

  /*!
   * Setup the \a PairCreator
   */
  virtual void setup();

  /*!
   * Create the pair linked list
   */
  virtual void createDistances();

  /*!
   * Clear all pair information. Particle positions have been updated.
   */
  virtual void invalidatePositions();

  /*!
   * Needed to initialise the old displacement stored in the tag
   */
  virtual void setupAfterParticleCreation();

  /*!
   * Return the maximal cutoff used for interactions
   * This means, this \a PairCreator subtracts the \a m_skinSize 
   */
  virtual double interactionCutoff() const {
    return ((Simulation*) ((Phase*) m_parent)->parent())->maxCutoff - m_skin_size;
  }


};


#endif
