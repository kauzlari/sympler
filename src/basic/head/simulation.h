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



#ifndef __SIMULATION_H
#define __SIMULATION_H 

#include <iostream>

#include "random.h"

#include "gen_f.h"
#include "meter.h"
#include "general.h"
#include "callable.h"
#include "controller.h"
#include "node_many_children.h"
#include "weighting_function.h"

using namespace std;


//---- Simulation ----

/*!
 * The main class that ties the simulation together
 */
class Simulation: public NodeManyChildren
{
protected:
  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * \a RandomNumberGenerator accessible by other modules and controlled by \a m_randomize
   */
  RandomNumberGenerator m_rng; 

  /*!
   * Create instances of the appropriate objects
   */
  virtual Node* instantiateChild(const string &name);

  /*!
   * Write restart file at the end of the simulation?
   */
  bool m_input_from_results;

  /*!
   * Name of this simulation
   */
  string m_name;

  /*!
   * Do we randomize \a m_rng ?
   */
  bool m_randomize;
  
  /*!
  * limt for number of loops to check for an endless loop during 
  * determination of \a Symbol stages
  */
  size_t m_stageSteps;

  /*!
   * Pointer to the \a Controller object
   */
  Controller *m_controller;

  /*!
   * Pointer to the \a Phase object
   */
  Phase *m_phase;

  /*!
   * List of \a WeightingFunction s
   */
  vector<WeightingFunction*> m_weighting_functions;

  /*!
   * List of \a Symbol s
   */
  vector<Symbol*> m_symbols;

  /*!
   * List of \a Callables s
   */
  vector<Callable*> m_callables;

  /*!
   * List of forces
   */
  vector<GenF*> m_forces;

  /*!
   * List of pair forces
   */
//   vector<GenF*> m_pair_forces;

  /*!
   * List of particle forces
   */
  vector<vector<GenF*>*> m_particle_forces;

  /*!
   * List of all other forces
   */
  vector<GenF*> m_other_forces;

  /*!
   * List of \a Meter s
   */
  vector<Meter*> m_meters;
		
public:
  /*!
   * Constructor
   */
  Simulation();

  /*!
   * Destructor
   */
  virtual ~Simulation();

  /*!
   * Run the simulation
   */
  void run();


  /*!
   * Returns a pointer to the \a Controller object
   */
  Controller* controller() {
    return m_controller;
  }


  /*!
   * Returns a pointer to the \a Phase object
   */
  Phase* phase() {
    return m_phase;
  }


  /*!
   * Returns the list of weighting functions
   */
  vector<WeightingFunction*> *weightingFunctions() {
    return &m_weighting_functions;
  }

  /*!
   * Returns the list of callables
   */
  vector<Callable*> *callables() {
    return &m_callables;
  }

  /*!
   * Returns the list of particle forces (for the parallel version)
   */
  vector<vector<GenF*>*> *particleForces() {
    return &m_particle_forces;
  }

  /*!
   * Returns the list of all other forces
   */
  vector<GenF*> *otherForces() {
    return &m_other_forces;
  }

  /*!
   * Returns the list of \a Meter s 
   */
  vector<Meter*> *meters() {
    return &m_meters;
  }

  /*!
   * Returns the list of \a Symbol  
   */
  vector<Symbol*> *symbols() {
    return &m_symbols;
  }

  /*!
   * Find the weighting function with name \a name
   * @param name Name of the weighting function
   */
  WeightingFunction *findWeightingFunction(const string &name);

  /*!
   * Find the force function with name \a name
   * @param name Name of the force function
   */
  GenF *findForceFunction(const string &name);

  GenF *searchForceWithClassName(const string &name, const string& className);


  /*!
   * Return the name of the simulation
   */
  const string &name() const {
    return m_name;
  }

  /*!
   * Pass program arguments to the \a Simulation
   */
  void readWithArg(const int& argc, char* argv[]);

  /*!
   * Write restart file
   */
  void writeWhenDone();

  /*!
   * Returns \a m_randomize 
   */
  bool randomize() {
    return m_randomize;
  }
 
  /*!
   * Different setups
   */
  virtual void setup();

  /*!
   * Maximum cut-off.
   */
  double maxCutoff;
	
  /*!
   * Call \a aboutToStart() for all \a Meter s
   */
  virtual void setupMeters();

  /*!
  * Setup the stages of the used symbols. NOTE: It is NOT enough to set the stages 
  * of the \a Symbol s in \a m_symbols.
  */
  virtual void setSymbolStages();

  /*!
  * Sort the stages of the used symbols. 
  */
  virtual void sortSymbolStages();


#ifdef _OPENMP
  /*! Set the offset to the copy vectors for the parallel version
   * Set the proper places to write in the copy vectors for each ValCalculator
   * Make a function out of it?
   */
  virtual void setupCopyVectors();
#endif
    

};

#endif
