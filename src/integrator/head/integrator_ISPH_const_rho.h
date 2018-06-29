/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#ifndef __INTEGRATOR_ISPH_CONST_RHO_H
#define __INTEGRATOR_ISPH_CONST_RHO_H

#include "integrator_position.h"

#include "function_pair.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorIISPHconstRho ----

/*!
 * Incompressible SPH \a Integrator based on the algorithm proposed by
 * S. Shao, E.Y.M. Lo / Advances in Water Resources 26 (2003) 787–800
 * A pressure Poisson equation (PPE) is SPH-siscretised and solved by 
 * relaxed Jacobi iteration.
 * To date (2017-10-12) it is still to be decided/tested, whether modified 
 * right-hand-sides (SPH-divergence of velocity instead of density finite 
 * difference in time) such as in 
 * B. Ataie-Ashtiani and G. Shobeyri, Int. J. Numer. Meth. Fluids 2008; 56:209–232
 * are better or not. Most likely the inplementation will allow to test both.
 */


class IntegratorISPHconstRho: public IntegratorPosition
{
protected:

  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * The mathematical expression aij' for the LHS of the SPH-
   * discretised pressure Poisson equation. See the documentation in 
   * \a setup() for further information.
   */
  string m_PPEexpr;
  
  /*!
   * The \a FunctionPair computing the user defined \a m_PPEexpr for
   * fluid-fluid particle pairs.
   */
  FunctionPair m_PPEfuncFluidFluid;
  
  /*!
   * The \a FunctionPair computing the user defined \a m_PPEexpr for
   * fluid-edge particle pairs.
   */
  FunctionPair m_PPEfuncFluidEdge;
  
  /*!
   * The \a FunctionPair computing the user defined \a m_PPEexpr for
   * fluid-wall particle pairs.
   */
  FunctionPair m_PPEfuncFluidWall;
  
  /*!
   * The \a FunctionPair computing the user defined \a m_PPEexpr for
   * edge-edge particle pairs.
   */
  FunctionPair m_PPEfuncEdgeEdge;
  
  /*!
   * The \a FunctionPair computing the user defined \a m_PPEexpr for
   * edge-wall particle pairs.
   */
  FunctionPair m_PPEfuncEdgeWall;
  
  /*!
   * Helper method for computation of a preliminary displacement increment due to the newest 
   * iterated pressure values
   */
  /* void computePairPressureIncrement(); */

  /*!
   * String listing all the other species which should be considered in the incompressibility algorithm
   */
  /* string m_speciesList; */

  /*!
   * helper list of all \a ColourPair s containing the incompressible species once
   */
  /* vector<ColourPair*> m_mixedColourPairs; */
  
  /*!
   * Name of the variable in the \a Particle tag for the generalised particle masses of all species
   */
  /* string m_genMassName; */

  /*!
   * Memory offset in the \a Particle tag for the generalised particle masses of all species
   */
  /* vector<size_t> m_genMassOffset; */

  /*!
   * The pointer to the implementation of the weighting function.
   */
  WeightingFunction* m_kernel;
  
  /*!
   * Name for the species of the edge particles placed directly on wall
   * boundaries.
   */
  string m_edgeSpecies;
  
  /*!
   * Name for the species of the wall particles placed in a range of 
   * 1 cutoff behind wall boundaries.
   */
  string m_wallSpecies;
  
  /*!
   * The name of the weighting function.
   */
  string m_kernelName;

  /*!
   * Name of the attribute in the \a Particle tag for the density
   */
  string m_densityName;

  /*!
   * Name of the list of bonded pairs that this Integrator 
   * automatically constructs for computation of interaction terms 
   * between two edge particles.
   */
  string m_edgeEdgeListName;

  /*!
   * Name of the list of bonded pairs that this Integrator 
   * automatically constructs for computation of interaction terms 
   * between edge particles and wall particles.
   */
  string m_edgeWallListName;

  /*!
   * Memory offset in the \a Particle tag for the density
   */
  size_t m_densityOffset;

  /*!
   * Colour of the species of the edge particles directly on 
   * the wall-boundaries
   */
  size_t m_edgeColour;
  
  /*!
   * Colour of the species of the particles behind the wall-boundaries
   */
  size_t m_wallColour;
  
  /*!
   * Memory offset in the \a Particle tag for the aii-matrix entries 
   * for the fluid (index 0) and edge particles (index 1)
   */
  size_t m_aiiOffset[2];

  /*!
   * Memory offset in the \a Particle tag for the preliminary advected 
   * density of the fluid (index 0), edge (index 1), and wall (index 2) 
   * species
   */
  size_t m_advDensityOffset[3];

  /*!
   * Memory offset in the \a Particle tag holding the old pressure 
   * within an iteration step for the fluid (index 0) and edge 
   * particles (index 1)
   */
  size_t m_pressureIterOldOffset[2];
  
  /*!
   * Memory offset in the \a Particle tag holding the new pressure 
   * within an iteration step for the fluid (index 0) and edge 
   * particles (index 1)
   */
  size_t m_pressureIterNewOffset[2];
  
  /*!
   * Memory offset in the \a Particle tag holding the pressure of the 
   * wall particles within an iteration step 
   */
  size_t m_pressureIterWallOffset;
  
  /*!
   * Index to access the list holding the bonded pairs of edge 
   * particles in the corresponding edge-edge \a ColourPair 
   */
  size_t m_edgeEdgeListIndex;
  
  /*!
   * Index to access the list holding the bonded pairs of edge and wall 
   * particles in the corresponding edge-wall \a ColourPair 
   */
  size_t m_edgeWallListIndex;
  
  /*!
   * Name of the attribute in the \a Particle tag of each species for the advected velocity
   * NOTE: See also the implemented help text in method init() 
   */
  /* string m_vAdvName; */

  /*!
   * Memory offset in the \a Particle tag for the preliminary advected velocity of all species
   */
  /* vector<size_t> m_vAdvOffset; */

  /*!
   * Memory offset in the \a Particle tag for the pair-contribution to the displacement increment 
   * due to the newest pressure
   */
  /* size_t m_pforcePairIncrOffset; */

  /*!
   * Memory offset in the \a Particle tag for the total displacement increment due to the newest pressure 
   * force. Note that \a m_pforcePairIncrOffset only stores the pair-contribution
   */
  /* size_t m_pforceIncrOffset; */

  /*!
   * Memory offset in the \a Particle tag for the density based on the newest iterated pressure
   */
  /* size_t m_iterDensityOffset; */

  /*!
   * Memory offset in the \a Particle tag to the pair contribution to 
   * the PPE matrix aij (j!=i) from the newest iterated pressure
   */
  size_t m_pairIterContribOffset[2];
  
  /*!
   * Name of the attribute in the \a Particle tag for the final incompressibility-preserving pressure
   */
  string m_pressureName;

  /*!
   * Memory offset in the \a Particle tag to the precomputable part of
   * the advected density for all \a Particle species
   */
  size_t m_advDensityPrecompOffset[3];

  /*!
   * Memory offset in the \a Particle tag to the interpolation 
   * normalisation for wall particles, when summing exclusively over
   * edge particles.
   */
  size_t m_wallEdgeNormalisationOffset;

  /*!
   * Memory offset in the \a Particle tag to the interpolation 
   * normalisation for edge and fluid particles.
   */  
  size_t m_advNormalisationOffset[2];

  /*!
   * Memory offset in the \a Particle tag to the acceleration caused
   * by the new pressure
   */  
  size_t m_pressureAccelFluidOffset;
  
  /*!
   * Memory offset in the \a Particle tag for the final incompressibility-preserving pressure
   */
  size_t m_pressureOffset[3];

  /*!
   * User defined maximal maximum relative density error for determination of the fulfillment of 
   * the incompressibility condition
   */
  /* double m_maxDensityError; */

  /*!
   * User defined maximal average relative density error for determination of the fulfillment of 
   * the incompressibility condition
   */
  /* double m_avgDensityError; */

  /*!
   * User defined maximum number of iterations for pressure computation. If the maximum is exceeded, 
   * the program is aborted.
   */
  size_t m_lMax;

  /*!
   * Fluid-fluid \a ColourPair for neighbour list managment
   */
  ColourPair* m_cpFluidFluid;
  
  /*!
   * Fluid-edge \a ColourPair for neighbour list managment
   */
  ColourPair* m_cpFluidEdge;
  
  /*!
   * Fluid-wall \a ColourPair for neighbour list managment
   */
  ColourPair* m_cpFluidWall;
  
  /*!
   * Edge-edge \a ColourPair for neighbour list managment
   */
  ColourPair* m_cpEdgeEdge;
  
  /*!
   * Edge-Wall \a ColourPair for neighbour list managment
   */
  ColourPair* m_cpEdgeWall;
  
  /*!
   * Relaxation parameter for the relaxed Jacobi iteration
   */
  double m_omega;
  
  /*!
   * Reference density for incompressibility condition
   */
  double m_rho0;
  
  /*! Convergence limit for the maximum error of the linear system
   *  solver for the PPE.
   */
  double m_epsilonMax;

  /*! Convergence limit for the average error of the linear system
   *  solver for the PPE.
   */
  double m_epsilonAvg;

  /*! Convergence limit for the difference of the maximum error of the linear system
   *  solver for the PPE between two iteration steps.
   */
  /* double m_epsilonMaxDelta; */

  /*! Convergence limit for the difference of the average error of the linear system
   *  solver for the PPE between two iteration steps.
   */
  /* double m_epsilonAvgDelta; */

  /*! Convergence limit for the maximum relative difference in
   *  iterated pressure values between two iteration steps. The
   *  relative difference is computed as (Pnew-Pold)/Pold. If Pold=0,
   *  (i.e. abs(Pold)<Peps, cf. attribute 'Peps') then Pnew is taken
   *  in the denominator.
   */
  double m_epsilonMaxRelDeltaP;
  
  /*! Threshold for computed pressures. Pressures below this threshold
   *  are treated as zero within the error estimate for the
   *  convergence check of the iterative solution.
   */
  double m_epsilonPZero;
  
  /*!
   * Helper that is not really necessary, but we need it at least now 
   * (20171206) due to reasons discussed at the initialisation in 
   * \a setup().
   */
  double m_kernelSelfContrib;
  
  /*!
   * Memory offset in the \a Particle tag meant to be used flexibly "as needed" for precomputed 
   * values. Currently (2017-01-26) it is used for only one value. Increased usage might save an 
   * insignificant amount of memory and make the code more messy
   */
  /* size_t m_precomputeOffset; */

  /*!
   * Helper that will be set to true if \a m_edgeSpecies and 
   * \a m_wallSpecies are set accordingly by user attributes 
   * 'edgeSpecies' and 'wallSPecies' 
   */
  bool m_usingWalls;

  /*!
   * Helper that will be set to true after this \a Integrator has
   * finished some precomputations that have to be performed only once 
   */
  bool m_precomputationDone;
  
  /*!
   * Helper function for creation of required particle colours.
   * ASSUMPTION: In any case \a IntegratorConstRho creates colour 0
   * which is interpreted as the fluid species. Depending on user 
   * input, 0 or 2 additional colours are created with colours 1 and 2 
   * in the latter case. This implies that no species and corresponding
   * colours may be created before the fluid species. And all 
   * additional species will have colours > 2. These particles from 
   * these additional species will be ignored by this 
   * \a IntegratorISPHConstRho.
   * @param colourTester Colour counter that is incremented by this method
   * @param species Species name of the colour to be created
   * @param attribute The user-attribute corresponding to the species 
   * in question. FIXME: This argument introduces redundancy and is, 
   * e.g., not safe against changes of the attribute name. The general 
   * FIX would be to store the attribute name in a variable as well 
   */
  void conditionalCreateColour(size_t& colourTester, string species, string attribute);

  /*!
   * Helper function that computes the pair contribution to the 
   * advected density from non-bonded pairs.
   * @param cp \a ColourPair of the contributing particle pairs
   */
  void advDensityPairSum(ColourPair* cp);

  /*!
   * Helper function that computes the pair contribution to the 
   * advected density from bonded pairs.
   * @param cp \a ColourPair of the contributing particle pairs
   * @param listIndex Index for access to the list of bonded pairs for 
   * the given \a ColourPair
   */
  void advDensPrecompConnectedPairContrib(ColourPair* cp, size_t listIndex);	

  /*!
   * Helper function that interpolates the "1" to obtain a 
   * normalisation for interpolated quantities.
   * @param cp \a ColourPair of the contributing particle pairs
   */
  void interpolateOne(ColourPair* cp);

  /*!
   * Helper function that computes the contributions of particle pairs 
   * to the diagonal matrix entries aii in the pressure Poisson 
   * equation.
   * @param cp \a ColourPair of the contributing particle pairs
   * @param PPEfunc \a FunctionPair computing the PPE pair expression 
   * defined in \a m_PPEexpr for the given \a ColourPair 
   */
  void aiiPairContrib(ColourPair* cp, const FunctionPair& PPEfunc);

  
  /*!
   * Helper function that initialises the pressure at the start of the 
   * iterative solving.
   * @param colour \a Colour of the particles to be initialised
   */
  void initPressure(size_t colour);

  /*!
   * Helper function that computes the pressure for wall particles from 
   * the pressure of edge the particles (assuming Neumann=0 BC).
   */
  void computeWallPressure();

  /*!
   * Helper function that computes the pair contribution aij, j!=i for 
   * the LHS of the PPE, based on the currently "new" pressure within 
   * the iteration.
   * @param cp The \a ColourPair of the species pair for which the 
   * contribution should be computed.
   * @param PPEfunc \a FunctionPair computing the PPE pair expression 
   * defined in \a m_PPEexpr for the given \a ColourPair 
   */
  void pairIterPContrib(ColourPair* cp, const FunctionPair& PPEfunc);

  
  /*!
   * Helper function that computes the total pair contribution 
   * aijPj i!=j to the pressure. The function checks itself if only for 
   * fluid or also for edge particles.
   */
  void totalPairContrib();

  /*!
   * Helper function that computes the final new pressure for the given
   * iteration step:
   * - Adds the RHS,  
   * - Multiplies with \a m_omega / aii
   * - Finally adds (1-omega)*Pold.
   * ASSUMPTION: The operations are performed on the data in the 
   * \a Particle tag accessed by \a m_pressureIterOldOffset for 
   * \a colour.
   * @param colour Colour of the particles this function should operate 
   * on
   */
  void newPressureIter(size_t colour);

  /*!
   * Performs operations required to compute the RHS of the PPE. The 
   * default behaviour in this abstract parent class is an empty 
   * function
   */
  virtual void computeRHS() {}
  
  /*!
   * Adds the RHS of the PPE to the storage of the new pressure in the
   * iteration. This is a pure virtual function in this class.
   * @param colour Colour of the particles this function should operate 
   * on
   */
  virtual void addRHStoNewPressure(size_t colour) = 0;
  
  /*!
   * Helper function that computes the total normalised l2-norm and the 
   * infinity norm for the residuals in the PPE solution   
   */
  void totalResiduals(double& maxRes, double& l2Res);

  /*!
   * Helper function that computes the maximum and squared l2-norm for 
   * residuals R_i`=Ri*dt^2 with 
   * R_i=RHS_i-AijPj = RHS_i - (AijPjpairs + AiiPi) 
   * for one colour.
   * @param colour Colour of the currently considered species of particles
   * @param maxRes Variable storing the largest residual found
   * @param sqL2Res Variable storing the squared l2-norm  
   */
  void residualsPerColour(size_t colour, double& maxRes, double& sqL2Res);

  /*!
   * Returns the largest relative pressure difference for the given 
   * iteration step l (Pnew) and the previous iteration step l-1 
   * (Pold) considering all relevant species. This function only 
   * collects and returns the maximum but does not any computation on 
   * its own
   */
  double returnMaxRelDeltaP();

  /*!
   * Computes the largest relative pressure difference for the given
   * iteration step l (Pnew) and the previous iteration step l-1 
   * (Pold) considering for the species given by \a colour. 
   * The formula is (Pnew-Pold)/Pold. If Pold<\a m_epsilonPZero
   * then Pnew is taken for the denominator.
   * @param colour Colour of the currently considered species of particles
   * @param maxRelDeltaP Variable holding the final result
   */
  void maxRelDeltaPPerColour(size_t colour, double& maxRelDeltaP);  
  
  /*!
   * Helper function for pressure force acceleration contribution FP/m 
   * of particles from other colours to the fluid particles.
   * We use the SPH-discretisation 
   * dv/dt = -nablaP/rho = -sum_j[mj(Pi/rhoi^2+Pj/rhoj^2)]nablaWij
   */
  void pressureForceIncrementFluidOther(ColourPair* cp);

  /*!
   * Helper function for pressure force acceleration contribution 
   * dt*FP/m of fluid particles to the fluid particles.
   * We use the SPH-discretisation 
   * dv/dt = -nablaP/rho = -sum_j[mj(Pi/rhoi^2+Pj/rhoj^2)]nablaWij
   */
  void pressureForceIncrementFluidFluid();

  /*!
   * Helper function for some variable initialisations. Should only be
   * called if \a m_usingWalls == true
   */
  void initialisationsWithWalls();

  /*!
   * Helper function for contributions to the advected density by wall 
   * and edge particles. Should only be called if 
   * \a m_usingWalls == true
   */
  void advDensityContribsWithWalls();


  void normalisationDenominatorsWithWalls();


  /*!
   * Helper function applying interpolation normalisations to the 
   * advected density for edge particles. Will only make changes if 
   * \a m_usingWalls == true
   */
  void normaliseAdvDensityWithWalls();

  
  /*!
   * Helper function for contributions to the diagonal part aii of the 
   * PPE by wall and edge particles. Should only be called if 
   * \a m_usingWalls == true
   */
  void aiiContribWithWalls();


  /*!
   * Helper function that computes the initial pressure for wall and 
   * edge particles. Should only be called if \a m_usingWalls == true
   */
  void initialPressureWithWalls();

    
  /*!
   * Helper function that initialises the wall and edge particle 
   * variables which are relevant for the pressure iteration. Should 
   * only be called if \a m_usingWalls == true
   */
  void initPressureIterWithWalls();


  /*!
   * Helper function that sets the final pressure to the last result 
   * obtained from the Jacobi iteration for wall and edge particles. 
   * Should only be called if \a m_usingWalls == true
   */
  void finalisePressureWithWalls();

  /*!
   * Helper function that computes the pair contribution within the 
   * Jacobi iteration from the edge and wall particles. 
   * Should only be called if \a m_usingWalls == true
   */
  void pairIterPContribWithWalls();

  
public:

  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
  IntegratorISPHconstRho(Controller *controller);

  /*!
   * Destructor
   */
  virtual ~IntegratorISPHconstRho();

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
  virtual void hitPos(const double& dt, const Particle* p, point_t &hit_pos,
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

  static const point_t s_dummyNullPoint;
  
};


#endif


