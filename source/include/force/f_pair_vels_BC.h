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


#ifndef __F_PAIR_VELS_BC_H
#define __F_PAIR_VELS_BC_H 

#include "f_pair_arbitrary_BC.h"
#include "pairdist.h"
#include "manager_cell.h"
// #include "function_pair.h"
#include "weighting_function.h"

using namespace std;

/*!
 * This is a completely general pair force FPV on a velocity v_i designed
 * for boundary interactions including Dirichlet boundary conditions
 * such that:
 * dv_i = FPV*dt
 *      = particleFactor_i(BV)*Sum_j(pairFactor_ij(BV)*weight_ij)*dt
 * where pairFactor_ij(BV) includes all pair contributions of the pair ij,
 * weight_ij represents the derivative of the used weighting function.
 * particleFactor_i(j)(BV) represents factors specific to particle i(j),
 * BV is the symbol for the boundary value (usually of some boundary-particles),
 * which may be used in the expressions as indicated.
 */
class FPairVelsBC : public FPairArbitraryBC
{
  protected:

    /*!
    * User defined string of the species forming the wall
    */
    string m_wallSpecies;

    /*!
     * Internal helper that sets the colour of the \a Wall \a Particle s according to \a m_wallSpecies
     */
    int m_wallColour;

    /*!
    * Internal helper that is set according to \a m_wallSpecies. It determines whether the \a Wall
    * \a Particle s are first or second in the \a Pairdist
    */
    bool m_wallIsSecond;

    /*!
     * Internal helper storing the \a Wall s in the non-periodic neighbouhood of each
     * \a Boundary \a Particle and the distance to it
     */
    vector<vector<pair<Wall*, double> >* >m_wallTable;

    /*!
     * Internal helper storing the \a Wall s in the periodic neighbouhood of each
     * \a Boundary \a Particle and the distance to it
     */
    vector<vector<pair<Wall*, double> >* >m_periodicWallTable;

/*!
 * Initialise the property list
 */
  void init();

  public:
 /*!
   * Constructor
   * @param simulation The \a Simulation object the force belongs to
  */
    FPairVelsBC(Simulation *simulation);

  /*!
     * Destructor
   */
    virtual ~FPairVelsBC();

  /*!
     * Setup this force, mainly the slots in memory
   */
    virtual void setup();

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
    * For each wall-particle, find the walls in range
    */
    virtual void setupAfterParticleCreation();
};
#endif
