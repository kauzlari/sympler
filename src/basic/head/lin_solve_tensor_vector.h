/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2016, 
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


#ifndef __LIN_SOLVE_TENSOR_VECTOR_H
#define __LIN_SOLVE_TENSOR_VECTOR_H

#include <fstream>

#include "general.h"
#include "function_pair.h"
#include "callable.h"

using namespace std;


/*!
 * Linear equation system solver. 
 * User defined input: tensor pair symbol for off-diagonal terms of matrix K, vector pair symbol holding solution v, vector pair symbol for right-hand side F. System to be solved for v is K.v=F. Diagonal term in a row of K is the negative sum of the off-diagonal terms of this row. 
 * Optional output: matrix K in ASCII or BINARY format.
 * FIXME: 1) The module does not yet solve but just outputs the matrix
 * FIXME: 2) Matrix is currently stored non-sparse
 */

class LinSolveTensorVector : public Callable
{
  protected:

  /*!
   * name of species 1.
   * FIXME: Flexible support of multiple species (and hence \a ColourPair s) missing. 
   * This Callable seems the easier way to realise this in the future than having a 
   * module for each \a ColourPair and finding a common place for storing the complete 
   * matrix.
   * Currently, two species can be defined. If they are diferent, the code should try
   * to create a matrix over the three possible species pairs. But even this is not 
   * implemented yet. Only two identical species are supported.
   */
  string m_species1;

  /*!
   * name of species 2.
   */
  string m_species2;

  /*!
   * Colour corresponding to \a m_species1
   */
  size_t m_colour1;

  /*!
   * Colour corresponding to \a m_species2
   */
  size_t m_colour2;

  /*!
   * name of the symbol representing the solution of the linear system
   */
  string m_solSymbolName;

  /*!
   * offset to the solution symbol with name \a m_solSymbolName for species1
   */
  size_t m_solSymbolOffset1;

  /*!
   * offset to the solution symbol with name \a m_solSymbolName for species2
   */
  size_t m_solSymbolOffset2;

#if 0
  /*! currently (2016-05-07, see directly below) realised by \a FunctionPair s*/
  /*!
   * name of the pair symbol representing the off-diagonal matrix entries 
   */
  string m_pairSymbolName;

  /*!
   * offset to the pair symbol with name \a m_pairSymbolName
   */
  size_t m_pairSymbolOffset;
#endif
  
  /*!
   * Cut-off radius for the pair summation
   */
  double m_cutoff;
  
  /*!
   * The mathematical pair-expression to be computed
   */
  string m_expression;
  
  /*!
   * The mathematical expression for \a m_1stparticleFactor11, \a m_1stparticleFactor12, \a m_1stparticleFactor22 
   */
  string m_1stPExpression;
  
  /*!
   * The mathematical expression for \a m_2ndparticleFactor11, \a m_2ndparticleFactor12, \a m_2ndparticleFactor22
   */
  string m_2ndPExpression;

  /*!
   * List of \a ColourPair s to work on. Has either 1 or 3 relevant entries
   */
  vector<ColourPair*> m_cps;

  /*!
   * List of \a FunctionPair s, each computing the user defined pair-expression for one \a ColourPair. Has either 1 or 3 relevant entries.
   */
  vector<FunctionPair*> m_functions;
  
  /*!
   * List of additional factors for the 1st particle for one \a ColourPair. Has either 1 or 3 relevant entries.
   */
  vector<FunctionPair*> m_1stParticleFactors;
  
  /*!
   * List of additional factors for the 2nd particle for one \a ColourPair. Has either 1 or 3 relevant entries.
   */
  vector<FunctionPair*> m_2ndParticleFactors;
  
  /*!
   * name of the symbol holding the right hand side of the equation
   */
  string m_RHSsymbolName;

  /*!
   * offset to the symbol with name \a m_RHSsymbolName for species1
   */
  size_t m_RHSsymbolOffset1;

  /*!
   * offset to the symbol with name \a m_RHSsymbolName for species2
   */
  size_t m_RHSsymbolOffset2;

  /*!
   * factor of 1 or -1 that makes the matrix symmetric or antisymmetric
   */
  int m_symmetry;

  /*!
   * how often to write the matrix to a file. 0 means never, 1 in every time step
   */
  size_t m_writeSystemEvery;

  /*!
   * If set to false, ascii files are produced, if set to true, binary files. In binary mode, of course no header is written
   */
  bool m_binary;

  /*!
   * name of file for writing the matrix.
   */
  string m_outputFileName;

  /*!
   * The output stream for the file with name \a m_outputFileName
   */
  ofstream m_s;

  /*!
   * The storage for the matrix
   */
  vector<double> m_mat;

  /*!
   * Helper for deciding when to write \a m_mat to a file
   */
  size_t m_countdown;

  /*!
   * Helper for numbering the output files
   */
  size_t m_step_counter;

  ///////////////////// Methods /////////////////////////

  /*!
   * Initialise the property list
   */
  virtual void init();
  
  /*!
   * writes header to ASCII output file
   */
  virtual void writeHeader();

  public:

    /*!
   * Constructor for the \a Node hierarchy
     */
    LinSolveTensorVector(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~LinSolveTensorVector();

  /*!
   * Starts the computations of this \a Callable
   */
  virtual void call(size_t timestep);

  /*!
         * Setup this \a Callable
   */
        virtual void setup();

  /*!
   * Additional setup
   */
  virtual void setupAfterParticleCreation();

};

#endif
