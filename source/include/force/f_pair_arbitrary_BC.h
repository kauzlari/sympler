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


#ifndef __F_PAIR_ARBITRARY_BC_H
#define __F_PAIR_ARBITRARY_BC_H 

#include "f_pair_wf.h"
#include "pairdist.h"
#include "manager_cell.h"
// #include "function_pair.h"
#include "weighting_function.h"
#include "function_pair_fixed.h"



using namespace std;

/*!
 * This is a completely general pair force FPV on a symbol x_i designed 
 * for boundary interactions including Dirichlet boundary conditions 
 * such that:
 * dx_i = FPV*dt
 *      = particleFactor_i(BV)*Sum_j(pairFactor_ij(BV)*weight_ij)*dt
 * where pairFactor_ij(BV) includes all pair contributions of the pair ij,
 * weight_ij represents the derivative of the used weighting function.
 * particleFactor_i(j)(BV) represents factors specific to particle i(j),
 * BV is the symbol for the boundary value (usually of some boundary-particles), 
 * which may be used in the expressions as indicated.
 */
class FPairArbitraryBC : public FPairWF
{
  protected:
  
  /*!
   * pair contribution to the force. It knows the symbol "BV" for the boundary value
   */
    FunctionPairFixed m_pairFactor;

   /*!
     * The mathematical expression for \a m_pairFactor
    */
    string m_pairFactorStr;
    
      /*!
     * An additional factor for the 1st particle. It knows the symbol "BV" for the boundary value
       */
    FunctionPairFixed m_1stparticleFactor;

    /*!
     * An additional factor for the 2nd particle. It knows the symbol "BV" for the boundary value
     */
    FunctionPairFixed m_2ndparticleFactor;

   /*!
     * The mathematical expression for \a m_1stparticleFactor
    */
    string m_1stPExpression;
    
   /*!
     * The mathematical expression for \a m_2ndparticleFactor
    */
    string m_2ndPExpression;

  /*!
     * Symmetry of the force
   */
    int m_symmetry;
  
    /*!
     * Initialise the property list
     */
    void init();

  public:
    /*!
   * Constructor
   * @param simulation The \a Simulation object the force belongs to
     */
    FPairArbitraryBC(Simulation *simulation);
    
    /*!
     * Destructor
     */
    virtual ~FPairArbitraryBC();

    /*!
     * Setup this force, mainly the slots in memory
     */
    virtual void setup();
};

#endif
