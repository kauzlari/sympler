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



#ifndef __TRIANGLE_CGMD_INTERPOLATION_H
#define __TRIANGLE_CGMD_INTERPOLATION_H

#include "general.h"
#include "callable.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

// the following definitions are done because I do not know how to turn range-checking 
// off for the GSL
#ifndef GSL_VECTOR_SET
#define GSL_VECTOR_SET(v, i, x) (v)->data[i*(v)->stride] = x
#endif
#ifndef GSL_VECTOR_GET
#define GSL_VECTOR_GET(v, i) (v)->data[i*(v)->stride]
#endif
#ifndef GSL_MATRIX_SET
#define GSL_MATRIX_SET(m, i, j, x) (m)->data[i * (m)->tda + j] = x
#endif
#ifndef GSL_MATRIX_GET
#define GSL_MATRIX_GET(m, i, j) (m)->data[i * (m)->tda + j]
#endif

using namespace std;

class Phase;
class Simulation;

/*!
 * A \a TriangleCGMDInterpolation interpolates atomic displacements, velocities, and forces on the nodes of a user-defined triangular mesh. Linear CGMD interpolation is used as described in Rudd&Broughton, Phys. Rev. B 72 144104
 * ASSUMPTION1: all original node-positions are INSIDE the box
 * ASSUMPTION2: the element sides are always smaller than box/2 in any direction
 * ASSUMPTION3: all elements lie in (x,y)-plane (is partially used, e.g., for 
 * shape functions and element area)
 */
class TriangleCGMDInterpolation: public Callable {

 protected:

  /*!
   * Time step when to activate this module
   */
  int m_activate_at;

  /*!
   * Time step when to deactivate this module
   */
  int m_deactivate_at;

  /*!
   * The interval when to call the module
   */
  size_t m_interval;

  /*!
   * The current time step
   */
  size_t m_step;

  /*!
   * Is the module active?
   */
  bool m_active;

  /*!
   * Species of the particles acting as nodes of the elements.
   */
  string m_nodeSpecies;

  /*!
   * Species of the particles delivering the data used for interpolation by the nodes.
   */
  string m_pSpecies;

  /*!
   * Colour corresponding to \a m_nodeSpecies
   */
  size_t m_nodeColour;

  /*!
   * Colour corresponding to \a m_pSpecies
   */
  size_t m_pColour;

  /*!
   * offset to the "source-particle's" symbol representing the 
   * slot of the corresponding "target-particle" ,i.e., "node"
   */
  size_t m_targetSlotOffset;

  /*!
   * name of the symbol representing the \a m_targetSlotOffset
   */
  string m_targetSlotName;

  /*!
   * offset to the particle's symbol representing the 
   * slot of the corresponding element it is in
   */
  size_t m_elementSlotOffset;

  //  /*!
  //   * name of the symbol representing the \a m_elementSlotOffset
  //   */
  //  string m_elementSlotName;

  /*!
   * Name of the file containing the nodal indices. It is assumed that 
   * these indices indicate slots of \a Particle s
   */
  string m_filename;

  /*!
   * Array storing the indices read from the file \a m_filename
   */
  vector<vector<size_t> > m_cornerIndices;

  /*!
   * Array storing the indices of the elements in \a m_cornerIndices for each node
   */
  vector<vector<int> > m_elementsOfNode;

  /*!
   * Array storing the three corners of each element in the shifted form without periodic images, i.e., the coordinates may differ from those of the nodes
   */
  vector<vector<point_t> > m_corners;

  /*!
   * Array storing the three corners of each element in their original form with (possibly) periodic images, i.e., the coordinates are the same as those of the nodes
   * FIXME: really necessary or node-coordinates sufficient?
   */
  vector<vector<point_t> > m_periodicCorners;

  /*!
   * The side-vectors of each element, pointing from one node to the other like this: 0->1->2->0
   * FIXME: do I really need to store this?
   */
  vector<vector<point_t> > m_sides;

  /*!
   * The vectors prependicular to the sides pointing INTO the element. This is necessary to detect whether a particle is inside an element 
   */
  vector<vector<point_t> > m_side_normals;

  /*!
   * The areas of each element
   */
  vector<double> m_areasTwice;

  /*
   * Array storing the products of summed contributions of particles i to the nodes m, n, i.e., N_mn = Sum_j(N_mj*N_nj)
   */
  vector<vector<double> > m_shapeMatrix;

  /*
   * Array storing in each time step the computed contributions of particles i to nodes m  according to the shape functions \a shapeF
   */
  vector<vector<double> > m_shapeF;


  // currently not needed
  //  /*
  //   * As \a m_shapeF, but for the initial configuration
  //   */
  //  vector<vector<double> > m_shapeF_0;

  /*
   * Array storing in each time step the computed contributions of particles i to nodes m  i.e., N_mi = Sum_n((Sum_j(N_mj*N_nj))^(-1))Nni, where the -1 indicates a matrix inversion
   */
  vector<vector<double> > m_f;

  // currently not needed
  //  /*
  //   * As \a m_f, but for the initial configuration
  //   */
  //  vector<vector<double> > m_f_0;

  // currently not needed
  //  /*
  //   * Array storing in each time step the difference \a m_f - \a m_f_0
  //   */
  //  vector<vector<double> > m_f_diff;

  /*!
   * Used for copying the matrix to be inverted with gsl_linalg_LU_decomp and gsl_linalg_LU_invert. 
   * This is OK/necessary because 
   * gsl_linalg_LU_decomp modifies it during LU decomposition
   */
  gsl_matrix* m_working_mat;
  
  /*!
   * Used for storing the inverse computed with gsl_linalg_LU_invert. 
   */
  gsl_matrix* m_inverse;
  
  /*!
   * The workspace for the computation of the LU decomposition with gsl_linalg_LU_decomp
   */
  gsl_permutation *m_permutation;

  /*!
   * Initial position of particles
   */
  vector<point_t> m_rInit;

  /*!
   * Name of memory offset to the interpolated displacement stored in the nodes
   */
  string m_displacementName;
  /*!
   * Name of memory offset to the non-locally interpolated momentum stored in the nodes
   */
  string m_momentum1Name;

  //  /*!
  //   * Name of memory offset to the locally interpolated momentum stored in the nodes
  //   */
  //  string m_momentum2Name;

  /*!
   * Name of memory offset to the non-locally interpolated force stored in the nodes
   */
  string m_force1Name;

  //  /*!
  //   * Name of memory offset to the locally interpolated force stored in the nodes
  //   */
  //  string m_force2Name;

  /*!
   * Memory offset to the interpolated displacement stored in the nodes
   */
  size_t m_displacementOffset;
  /*!
   * Memory offset to the non-locally interpolated momentum stored in the nodes
   */
  size_t m_momentum1Offset;

  //  /*!
  //   * Memory offset to the locally interpolated momentum stored in the nodes
  //   */
  //  size_t m_momentum2Offset;

  /*!
   * Memory offset to the non-locally interpolated force stored in the nodes
   */
  size_t m_force1Offset;

  //  /*!
  //   * Memory offset to the locally interpolated force stored in the nodes
  //   */
  //  size_t m_force2Offset;

  /*!
   * If set to true, the interpolation matrix will not be computed as a 
   * function of the initial (assumed as equilibrium) atomic positions 
   * but as a function of the actual positions at every time step. This 
   * is of course more time-consuming because it requires the 
   * corresponding matrix inversion at every time step and not just 
   * once at the beginning.
   */
  bool m_transient;


  /*!
   * Should the matrix of shape functions be written to a file?
   */
  bool m_writeNN;

  /*!
   * Internal counter for deciding when to write the matrix of shape functions
   */
  size_t m_NNstepCounter;

  /*!
   * Internal counter for the matrix of shape functions files 
   */
  size_t m_NNfileCounter;

  /*!
   * User input for deciding how often to write the matrix of shape functions
   */
  size_t m_writeNNevery;

  /*!
   * File name for storing the matrix of shape functions in binary format. A counter and the ending '.bin' are added automatically.
   */
  string m_NNname;

  /*!
   * Should the matrix of shape functions be written to a file?
   */
  bool m_write_f;

  /*!
   * Internal counter for deciding when to write the matrix of shape functions
   */
  size_t m_fStepCounter;

  /*!
   * Internal counter for the CGMD interpolation matrix 
   */
  size_t m_fFileCounter;

  /*!
   * User input for deciding how often to write the CGMD interpolation matrix
   */
  size_t m_write_fEvery;

  /*!
   * File name for storing the CGMD interpolation matrix in binary format. A counter and the ending '.bin' are added automatically.
   */
  string m_fName;

  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Helper: Should the module apply periodic boundary conditions?
   */
  bool m_periodic;

  /*!
   * Helper: for applying periodic boundary conditions
   */
  point_t m_boxSize;

  /*!
   * Helper: was the interpolation matrix of the initial configuration already stored?
   */
  bool m_zeroStored;
  
  /*!
   * Function for debugging: print m_f 
   */
  void debugPrintM_F(size_t node) const;

  /*!
   * Function for debugging: print the sum of the node rows m_f for each particle 
   */
  void  debugSumPerPart() const;

  /*!
   * Function for debugging: print the sum of the particle columns m_f for each node 
   */
  void  debugSumPerNode() const;

  /*!
   * Function for debugging: print the sum of all entries of m_f 
   */
  void  debugTotalSum() const;

  /*!
   * Computation of the linear shape function of local node n of the triangular 
   * element at position pr
   * ASSUMPTION: pr is inside the element!!!
   * ASSUMPTION3 is also used !!!
   */
  inline double shapeF(size_t el, size_t n, const point_t& pr) const {
    const vector<point_t>& elCorner = m_corners[el];
    if(n == 0)
    // node0: [x1y2-x2y1 + (y1-y2)x + (x2-x1)y]/2Area
      return 
	( 
	 elCorner[1][0]*elCorner[2][1] - elCorner[2][0]*elCorner[1][1] // x1y2-x2y1
	 +(elCorner[1][1]-elCorner[2][1])*pr.x // + (y1-y2)x
	 +(elCorner[2][0]-elCorner[1][0])*pr.y // + (x2-x1)y
	 ) / m_areasTwice[el];
    if(n==1)
      // node1: [x2y0-x0y2 + (y2-y0)x + (x0-x2)y]/2Area
      return
	( 
	 elCorner[2][0]*elCorner[0][1] - elCorner[0][0]*elCorner[2][1] // x2y0-x0y2
	 +(elCorner[2][1]-elCorner[0][1])*pr.x // + (y2-y0)x
	 +(elCorner[0][0]-elCorner[2][0])*pr.y // + (x0-x2)y
	 ) / m_areasTwice[el];
    if(n==2)
      // node2: [x0y1-x1y0 + (y0-y1)x + (x1-x0)y]/2Area
      return
	( 
	 elCorner[0][0]*elCorner[1][1] - elCorner[1][0]*elCorner[0][1] // x0y1-x1y0
	 +(elCorner[0][1]-elCorner[1][1])*pr.x // + (y0-y1)x
	 +(elCorner[1][0]-elCorner[0][0])*pr.y // + (x1-x0)y
	 ) / m_areasTwice[el];
    throw gError("TriangleCGMDInterpolation::shapeF"+FILE_INFO, "called with forbidden node number " + ObjToString(n) + "! Allowed are {0, 1, 2}!");
    return HUGE_VAL;
  }
  
  /*!
   * Computation of the linear shape function of node0 of the triangular 
   * element at position pr
   * ASSUMPTION3 is used !!!
   */
  inline double shapeF0(size_t el,  const point_t& pr) const {
    const vector<point_t>& elCorner = m_corners[el];

    // node0: [x1y2-x2y1 + (y1-y2)x + (x2-x1)y]/2Area
    return 
      (
       elCorner[1][0]*elCorner[2][1] - elCorner[2][0]*elCorner[1][1] // x1y2-x2y1
       +(elCorner[1][1]-elCorner[2][1])*pr.x // + (y1-y2)x
       +(elCorner[2][0]-elCorner[1][0])*pr.y // + (x2-x1)y
       ) / m_areasTwice[el];
  }

  /*!
   * Computation of the linear shape function of node1 of the triangular 
   * element at position pr
   * ASSUMPTION3 is used !!!
   */
  inline double shapeF1(size_t el,  const point_t& pr) const {
    const vector<point_t>& elCorner = m_corners[el];

    // node1: [x2y0-x0y2 + (y2-y0)x + (x0-x2)y]/2Area
    return
      ( 
       elCorner[2][0]*elCorner[0][1] - elCorner[0][0]*elCorner[2][1] // x2y0-x0y2
       +(elCorner[2][1]-elCorner[0][1])*pr.x // + (y2-y0)x
       +(elCorner[0][0]-elCorner[2][0])*pr.y // + (x0-x2)y
       )/ m_areasTwice[el];
  }

  /*!
   * Computation of the linear shape function of node 2 of the triangular 
   * element at position pr
   * ASSUMPTION3 is used !!!
   */
  inline double shapeF2(size_t el,  const point_t& pr) const {
    const vector<point_t>& elCorner = m_corners[el];

    // node2: [x0y1-x1y0 + (y0-y1)x + (x1-x0)y]/2Area
    return 
      (
       elCorner[0][0]*elCorner[1][1] - elCorner[1][0]*elCorner[0][1] // x0y1-x1y0
       +(elCorner[0][1]-elCorner[1][1])*pr.x // + (y0-y1)x
       +(elCorner[1][0]-elCorner[0][0])*pr.y // + (x1-x0)y
       )/ m_areasTwice[el];
  }

 public:
  /*!
   * Constructor
   * @param sim Pointer to the main simulation object
   */
  TriangleCGMDInterpolation(Simulation* sim);
  
  /*!
   * Destructor
   */
  virtual ~TriangleCGMDInterpolation();
  
  /*!
   * Starts the computations
   */
  virtual void call(size_t timestep);
  
  /*!
   * Does the actual computations
   */
  virtual void compute();
  
  /*!
   * Interpolates particle information to nodes. What is interpolated is currently (2011-03-21) hard-coded
   */
  virtual void interpolate() const;
  
  /*!
   * Setup this \a TriangleCGMDInterpolation
   */
  virtual void setup();
  
  /*!
   * Additional setup and execution of \a call(0)
   */
  virtual void setupAfterParticleCreation();
    
};

#endif
