/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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




#ifndef __INTEGRATOR_IISPH_CONST_RHO_H
#define __INTEGRATOR_IISPH_CONST_RHO_H

#include "integrator_position.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorIISPHconstRho ----

/*!
 * Implicit incompressible SPH \a Integrator based on the algorithm proposed by
 * Ihmsen at al. (IEEE Transactions on Visualization and Computer Graphics 20, 426 (2014))
 */

class IntegratorIISPHconstRho: public IntegratorPosition
{
protected:
  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Helper method for setup of memory offsets in the \a Particle tag for newly created attributes
   * The attributes used here are all either initialised by hand or must indeed be persistent. Hence
   * the default value below.
   */
  size_t addNewAttr(size_t colour, string symbolName, DataFormat::datatype_t datatype, bool persistency = true);

  /*!
   * Helper method for computation of a preliminary displacement increment due to the newest 
   * iterated pressure values
   */
  void computePairPressureIncrement();

  /*!
   * String listing all the other species which should be considered in the incompressibility algorithm
   */
  string m_speciesList;

  /*!
   * helper list of all \a ColourPair s containing the incompressible species once
   */
  vector<ColourPair*> m_mixedColourPairs;
  
  /*!
   * Name of the variable in the \a Particle tag for the generalised particle masses of all species
   */
  string m_genMassName;

  /*!
   * Memory offset in the \a Particle tag for the generalised particle masses of all species
   */
  vector<size_t> m_genMassOffset;

  /*!
   * The name of the weighting function.
   */
  string m_kernelName;

  /*!
   * The pointer to the implementation of the weighting function.
   */
  WeightingFunction* m_kernel;

  /*!
   * Name of the attribute in the \a Particle tag for the density
   */
  string m_densityName;

  /*!
   * Memory offset in the \a Particle tag for the density
   */
  size_t m_densityOffset;

  /*!
   * Memory offset in the \a Particle tag for the dii-matrix \a point_t entries of the incompressible species
   */
  size_t m_diiOffset;

  /*!
   * Memory offset in the \a Particle tag for the aii-matrix double entries of the incompressible species
   */
  size_t m_aiiOffset;

  /*!
   * Memory offset in the \a Particle tag for the preliminary advected density of the incompressible species
   */
  size_t m_advDensityOffset;

  /*!
   * Name of the attribute in the \a Particle tag of each species for the advected velocity
   * NOTE: See also the implemented help text in method init() 
   */
  string m_vAdvName;

  /*!
   * Memory offset in the \a Particle tag for the preliminary advected velocity of all species
   */
  vector<size_t> m_vAdvOffset;

  /*!
   * Memory offset in the \a Particle tag for the pair-contribution to the displacement increment 
   * due to the newest pressure
   */
  size_t m_pforcePairIncrOffset;

  /*!
   * Memory offset in the \a Particle tag for the total displacement increment due to the newest pressure 
   * force. Note that \a m_pforcePairIncrOffset only stores the pair-contribution
   */
  size_t m_pforceIncrOffset;

  /*!
   * Memory offset in the \a Particle tag for the density based on the newest iterated pressure
   */
  size_t m_iterDensityOffset;
  
  /*!
   * Memory offset in the \a Particle tag for the new newest pressure computed in an iteration
   * NOTE: The offset will be swapped with \a m_pressureIterOffsetOld in each iteration step
   */
  size_t m_pressureIterOffsetNew;

  /*!
   * Memory offset in the \a Particle tag for the new old pressure still needed in an iteration
   * NOTE: The offset will be swapped with \a m_pressureIterOffsetNew in each iteration step
   */
  size_t m_pressureIterOffsetOld;

  /*!
   * Name of the attribute in the \a Particle tag for the final incompressibility-preserving pressure
   */
  string m_pressureName;

  /*!
   * Memory offset in the \a Particle tag for the final incompressibility-preserving pressure
   */
  size_t m_pressureOffset;

  /*!
   * User defined maximal maximum relative density error for determination of the fulfillment of 
   * the incompressibility condition
   */
  double m_maxDensityError;

  /*!
   * User defined maximal average relative density error for determination of the fulfillment of 
   * the incompressibility condition
   */
  double m_avgDensityError;

  /*!
   * User defined maximum number of iterations for pressure computation. If the it is exceeded, 
   * the program is aborted.
   */
  size_t m_lMax;

  /*!
   * Relaxation parameter for the relaxed Jacobi iteration
   */
  double m_omega;
  
  /*!
   * Reference density for incompressibility condition
   */
  double m_rho0;
  
  /*!
   * Memory offset in the \a Particle tag meant to be used flexibly "as needed" for precomputed 
   * values. Currently (2017-01-26) it is used for only one value. Increased usage might save an 
   * insignificant amount of memory and make the code more messy
   */
  size_t m_precomputeOffset;
  
public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorIISPHconstRho(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorIISPHconstRho();

  /*!
   * Setup for this \a Integrator
   */
  virtual void setup();

  /*!
   * Initialize temporary fields and clear all forces
   */
  virtual void isAboutToStart();

  /*!
   * Position integration in collaboration with the \a Cell s
   */
  virtual void integrateStep1();

  /*!
   * Integration of the velocity by iterative pressure correction = IISPH alrorithm
   */
  virtual void integrateStep2();

  /*!
   * Integration of the position
   */
  virtual void integratePosition(Particle* p, Cell* cell);

  /*!
   * Prediction of the velocity; 
   * called within integrateStep1; does nothing for this \a IntegratorPosition
   */
  virtual void integrateVelocity(Particle* p) {}

  /*!
   * Solves the equation that checks for hits.
   * NOTE: currently (2017-01-05) not supported and hence throws exception since
   * \a integratePosition should not trigger this function any way 
   * Most kinds of 
   * collisons and subsequent reflections will violate incompressibility anyway. This 
   * means wall collisions are not detected and hence wall penetration and loss of 
   * particles may occur!
   */
  virtual void solveHitTimeEquation(WallTriangle* wallTriangle, const Particle* p, const point_t
				    &force, vector<double>* results);

  /*!
   * Returns the actual hit position at a \a WallTriangle in the argument \a hit_pos. 
   * The function is used in \a WallTriangle.
   * NOTE: For the same reasons as in \a solveHitTimeEquation, currently (2017-01-05), this 
   * function should never be called and hence throws an exception
   */
  virtual void hitPos(double dt, const Particle* p, point_t &hit_pos,
		      const point_t &force);

#ifdef _OPENMP
  /*!
   * Returns a characteristic string for the integrated degrees of freedom (DOFs)
   */
  virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

#endif

  static const point_t dummyNullPoint;
  
};


#endif


