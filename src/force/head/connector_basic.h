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



#ifndef __CONNECTOR_BASIC_H
#define __CONNECTOR_BASIC_H

#include "gen_connector.h"
#include "pairdist.h"
#include "simulation.h"
#include "manager_cell.h"
#include "function_pair.h"

using namespace std;

// class FunctionPair;

/*!
 * A basic connection between two particles, defining a pair force
 */
class ConnectBasic : public GenConnector
{
protected:

  /*!
  * The mathematical expression for the connector-force
  */
  FunctionPair m_pairFactor;

  /*!
   * Initialize the property list
   */
  void init();

public:
  /*!
   * Constructor for \a Node hierarchy
   */
  ConnectBasic(Simulation *simulation);

  /*!
  * Standard constructor will throw an exception
  */
  ConnectBasic();

  /*!
   * Destructor
   */
  virtual ~ConnectBasic();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  /*!
   * Compute the force
   * @param force_index The index for the memory slot to save the current force
     */
  virtual void computeForces(int force_index);

#ifndef _OPENMP
  /*!
   * Compute the force
   * @param force_index The index for the memory slot to save the current force
   * @param pair is the current Pairdist, the force acts on.
   */
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

  /*!
   * Compute the force
   * @param force_index The index for the memory slot to save the current force
   * @param part is the current Particle this force acts on.
   */
  virtual void computeForces(Particle* part, int force_index);

  /*!
   * Setup shortly before simulation starts
   */
  virtual void setup();
};

#endif
