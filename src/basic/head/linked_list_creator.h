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



#ifndef __LINKED_LIST_CREATOR_H
#define __LINKED_LIST_CREATOR_H

using namespace std;


#include <list>

#include "cell.h"
#include "gen_f.h"
#include "pairdist.h"
#include "pair_list.h"
#include "wall_container.h"
#include "pair_creator.h"
#include "node_many_children.h"


class LinkedListCreator: public PairCreator
{
protected:

  /*!
   * Initialize
   */
  void init();


public:

  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   */
  LinkedListCreator(Phase* p);

  /*!
   * Destructor
   */
  virtual ~LinkedListCreator();

  /*!
   * Setup the linked list creator
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

};


#endif
