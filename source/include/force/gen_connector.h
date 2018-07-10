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



#ifndef __G_CONNECTOR_H
#define __G_CONNECTOR_H

#include "gen_f.h"
#include <map>
#include <list>
#include <utility>
#include <fstream>

#include "misc.h"
#include "pairdist.h"
#include "pair_list.h"
#include "integrator_position.h"
#include "integrator_velocity_verlet.h"
#include "integrator_velocity_verlet_disp.h"
#include "integrator_static.h"
#include "colour_pair.h"

using namespace std;

class Simulation;

/*!
 * Base class for all connection forces.
 * It also handles
 * the creation of the pair list for which the connection is relevant.
 */
class GenConnector : public GenF
{
protected:

  /*!
   * The species of the two in the force participation particles.
   */
  pair<string, string> m_species;

  /*!
   * The colour pair belonging to the species combination.
   */
  ColourPair *m_cp;

  /*!
   * Inititialize the property list.
   */
  void init();
   
public:
  /*!
   * Constructor
   */
  GenConnector();	
  GenConnector(Simulation *simulation);
  GenConnector(Simulation *simulation,ColourPair* cp);

  /*!
   * The list of connections this connector is working on. It must be requested from the corresponding \a ColourPair linked in \a m_cp
   */
  PairList* m_connectedList;


  /*!
   * Destructor
   */
  virtual ~GenConnector();

  /*!
   * setup this class
   */
  virtual void setup();

  /*!
   * Extra setup when all particles have been created
   */
  virtual void setupAfterParticleCreation();


  /*!
   * Adds a particle pair to the connected list
   */
  // OLD: 2 lines; now in ColourPair
/*   void addPairToConnection(Particle *first_p, Particle *second_p); */
/*   void updateDistances(); */

};

#endif







