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



#ifndef __CONTROLLER_H
#define __CONTROLLER_H

#include <iostream>

#include "integrator.h"
#include "node_many_children.h"


using namespace std;



class Simulation;

/* Controller has integrators as its children.
   Minimum is one integrator (for position and momentum). */

/*!
 * A \a Controller handles the main simulation loop.
 * For each time step is calls:
 *   1. Print the current time step to screen
 *   2. Call \a Meter s
 *   3. Integrate, which means
 *        a) Call \a integrateStep1 of all \a Integrator s
 *        b) Clear all temporary particle/pair information
 *        c) Compute forces
 *        d) Call \a integrateStep2 of all \a Integrator s
 *   4. Call \a Callable s
 *
 * A \a Controller has \a Integrator s as its children.
 */
class Controller: public NodeManyChildren
{
protected:
  /*!
   * Total number of time steps
   */
  int m_timesteps;

  /*!
   * Simulation time step
   */
  double m_dt;

  /*!
   * Current simulation time
   */
  double m_t;

  /*!
   * How often to print out the simulation progress
   */
  int m_statusEvery;

  /*!
   * The current force index in the force history table (the Newtonian force only)
   */
  int m_force_index;

  /*!
  * A list where \a Node s may register for being called for an additional setup 
  * after particle creation
  */
  list<Node*> m_toSetupAfterParticleCreation;
  
  /*!
  * A list where \a Node s may register for being called for some precomputation operation before the s\a Symbol ("stage 1") and \a Force computations
  */
  list<Node*> m_precomputers;
  
  /*!
  * A list where \a Node s may register for being called for some precomputation operation before the s\a Symbol computations in "stage 0".
  */
  list<Node*> m_precomputers_0;
  
  /*!
   * Register additional properties with the property list
   */
  void init();

  /*!
   * Integrate the equation of motion, which means
   *        a) Call \a integrateStep1 of all \a Integrator s
   *        b) Clear all temporary particle/pair information
   *        c) Compute forces
   *        d) Call \a integrateStep2 of all \a Integrator s
   */
  void integrate();

  /*!
   * Look for integrators in the \a Integrator database and instantiate if
   * found.
   */
  virtual Node* instantiateChild(const string &name);

  void runSymbols();

  void runSymbols_0();

public:
static /*double*/int time_for_parallel;
static /*double*/int time_for_parallel1;
static /*double*/int time_for_parallel2;
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  Controller(Simulation *simulation);

  /*!
   * Destructor
   */
  ~Controller();

  /*!
   * \a setup initializes all the variables. This ensures object
   * creation of all entities in the XML file is complete when setup is being called.
   * It also invokes setup for all child objects.
   */
  virtual void setup();

  /*!
   * Find the integrator \a name for species \a species
   * @param name Name of the integrator
   * @param species Species this integrator acts on
   */
  Integrator *findIntegrator(const string &name, const string &species) {
    return (Integrator*) findChild(name + "(" + species + ")");
  }
	
  /*!
   * Run the simulation
   */
  void run();

  /*!
   * Return a pointer to the list of integrators. Fixme!!! If not obsolete,
   * at least this is not elegant
   */
  list<Node*>* integrators() {
    return &m_children;
  }

  /*!
   * Return the total number of time steps
   */
  int timesteps() const {
    return m_timesteps;
  }
  
  /*!
   * Return the current time
   */
  double time() const {
    return m_t;
  }

  /*!
   * Return the time step
   */
  double dt() const {
    return m_dt;
  }

  /*!
   * Return the force index for access to the force history table
   */
  int forceIndex() const {
    return m_force_index;
  }
  
  /*!
  * Adds the \a Node to a list for being called before \a Symbol ("stage 1") and \a Force computations
  */
  void registerForPrecomputation (Node* callable);
  
  /*!
  * Adds the \a Node to a list for being called before the \a Symbol computations of "stage 0"
  */
  void registerForPrecomputation_0 (Node* callable);
  
  /*!
  * Adds the argument to a list for being called after particle creation
  */
  void registerForSetupAfterParticleCreation (Node* callable);
/*   { */
/*     m_toSetupAfterParticleCreation.push_back(callable);  */
/*   } */


 private:

  double m_pairCreateTime;
  double  m_pairForceTime;
  double  m_otherForceTime;
  double  m_integrateTime;
  double m_initTime;
  double m_partSymbolsTime;
  double m_bondedSymbolsTime;
  double m_nonBondedSymbolsTime;
  double m_SpecialTime;

};
#endif
