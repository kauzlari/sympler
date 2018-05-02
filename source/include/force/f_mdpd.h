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



#ifndef __F_MDPD_H
#define __F_MDPD_H 

#include "f_pair_wf.h"
#include "pairdist.h"
#include "manager_cell.h"
#include "function_fixed.h"
#include "weighting_function.h"


using namespace std;


/*!
 * Implementation of general MDPD forces with influence from the wall.
 */

class FMDPD : public FPairWF
{
protected:
  /*!
   * For passing information to the threads.
   * The current index in the force history.
   */
  size_t m_force_index;

  /*!
   * Offset for access to the local density for each species
   */
  pair<size_t, size_t> m_density_offset;

  /*!
   * Name of the local density to be used
   */
  string m_rhoSymbol;

  /*!
   * The dpsi/dn function
   */
  FunctionFixed m_dpsi_dn;

  /*!
   * This is the maximum density which is used in calculations.
   * Larger ones will be reduced to it.
   */
  double m_density_cutoff;

  /*!
   * Will the local density be computed over all ColourPairs (CP) or only over the CP
   * corresponding to the chosen species?
   */
  bool m_oneProp;

  void init();

public:
  /*!
   * Constructor.
   */
  FMDPD(Simulation *simulation);

  /*!
   * Destructor.
   */
  virtual ~FMDPD();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  virtual void computeForces(int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif

  virtual void computeForces(Particle* part, int force_index);

  virtual void setup();
};

#endif
