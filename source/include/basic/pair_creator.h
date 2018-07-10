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




#ifndef __PAIR_CREATOR_H
#define __PAIR_CREATOR_H

using namespace std;


#include <list>

/* #include "cell.h" */
#include "node.h"
#include "gen_f.h"
#include "pairdist.h"
#include "pair_list.h"
#include "wall_container.h"
#include "phase.h"
#include "node.h"
#include "simulation.h"


class ParticleCreator;
class ColourPair;


class PairCreator: public Node
{
protected:

  /*!
   * Name of the pair creator
   */
  string m_pairCrName;

  /*!
   * Are the distances valid or do they need to be recalculated next time?
   * Fixme!!! quick'n'dirty
   */
  bool m_valid_dist;

  /*!
   * Initialise the property list
   */
  virtual void init();


public:

  /*!
   * Constructor
   * @param p Pointer to the parent phase object
   */
  PairCreator(Phase* p);

  /*!
   * Destructor
   */
  virtual ~PairCreator();

  /*!
   * At the moment doesn't really do anything
   */
  virtual void setup();

  /*!
   * Create the pair list
   */
  virtual void createDistances();

  /*!
   * Sets the moment when we have to recreate the pair list
   */
  virtual void invalidatePositions();

  /*!
   * Lets the \a PairCreator perform an action 
   * when a pair of particles has been created
   */
  virtual void pairAddedFree(Pairdist* pd);

  /*!
   * Lets the \a PairCreator perform an action 
   * when a (frozen) pair of particles has been created
   */
  virtual void pairAddedFrozen(Pairdist* pd);

  /*!
   * The counter is needed to tell the created new pairs in which Pairlist they belong.
   */
  static int counterTN;

  /*!
   * Return whether the distances between the particles need to be created
   */
  virtual bool& validDist() {
    return m_valid_dist;
  }

  /*!
   * Return the maximal cutoff used for interactions
   */
  virtual double interactionCutoff() const {
    return ((Simulation*) ((Phase*) m_parent)->parent())->maxCutoff;
  }

};

//---- Factories ----

class PairCreator_Factory: public SmartEnum<PairCreator_Factory>
{
public:
  virtual PairCreator *instantiate(Phase *phase) const = 0;

protected:
  PairCreator_Factory(const string &name)
      : SmartEnum<PairCreator_Factory>(name) { }
};


template <class T>
class PairCreator_Register: public PairCreator_Factory
{
public:
    PairCreator_Register(const string &name)
    : PairCreator_Factory(name) { }

    virtual PairCreator *instantiate(Phase *phase) const;
};



//---- Inline functions ----

template <class T>
inline PairCreator *PairCreator_Register<T>::instantiate(Phase *phase) const
{
    return new T(phase);
}


#endif

