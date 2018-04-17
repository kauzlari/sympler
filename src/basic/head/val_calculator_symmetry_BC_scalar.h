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


#ifndef __VAL_CALCULATOR_SYMMETRYBC_SCALAR_H
#define __VAL_CALCULATOR_SYMMETRYBC_SCALAR_H

#include "val_calculator_BC.h"

// #include "wall.h"

class Wall;

//---- newclass -----------------------------------------------------------

/*!
 * Saves the pair-specific value of a quantity of the boundary particle used
 * for applying a symmetry boundary condition (BC) in each pair of particles.
 * This class implements symmetry BCs for an arbitrary scalar.
 */
class ValCalculatorSymmetryBCScalar : public ValCalculatorBC
{
  protected:

  /*!
   * Name of the scalar the BC is computed for
   */
  string m_scalarName;

  /*!
   * Memory offset to the scalar descibed by \a m_scalarName
   */
  pair<size_t, size_t> m_scalarOffset;
    
  /*!
   * Name of the gradient used for computing the BC
   */
  string m_gradientName;

  /*!
   * Memory offset to the gradient descibed by \a m_gradientName
   */
  size_t m_gradientOffset;

  /*!
   * Helper store the right wall for the current \a Pairdist
   */
  Wall* m_currentWall; 
    
    /*!
     * Initialise the property list
     */
    virtual void init();
  
    /*!
     * Helper function for returning a copy of itself
     */
    virtual ValCalculatorPair* copyMySelf()
    {
      return new ValCalculatorSymmetryBCScalar(*this); 
    }
    
    /*!
     * Find the wall-segment where the rij-vector is hitting and the distance 
     * of i and j to the \a Wall
     * Additionally this function sets \a m_currentWall
     * @param pD \a Pairdist for which to compute the distances
     * @param innerDist variable for distance of inner particle to wall
     * @param outerDist variable for distance of wall particle to wall
     */ 
    virtual void computeDists(Pairdist* pair, double& innerDist, double& outerDist); //FIXME: inline ?

    /*!
     * Returns the strings of those \a Symbols that the given class depends on
     * due to hard-coded reasons (not due to runtime compiled expressions).
     * @param usedSymbols List to add the strings to.
     */
    virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
    {
      throw gError("ValCalculatorSymmetryBCScalar", "This exception was created while this module was inactive (not compiled), and it is meant to force you to at least consider the creation of a common parent class together with ValCalculatorDicrichletBCScalar");
      usedSymbols.push_back(m_scalarName);
    }

    
  public:
  /*!
   * Constructor
   * @param symbol The symbol name
   */
    ValCalculatorSymmetryBCScalar(string symbol);

  /*!
     * Constructor for Node hierarchy
   */
    ValCalculatorSymmetryBCScalar(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorSymmetryBCScalar() {
    }
        
    /*!
     * Compute the boundary value for the \a Pairdist \a pD
     * @param pD \a Pairdist for which to compute the boundary value
     */ 
    virtual void compute(Pairdist* pair); //FIXME: inline ?

  /*!
     * Returns "symmetryBCScalar"
   */
    virtual string myName() {
      return "symmetryBCScalar";
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();

    /*!
     * Find the right stage for the computation right at the beginning of a timestep.
     */
    bool findStage_0();

    /*!
     * Find the right stage for the computation before the forces.
     */
    bool findStage();

};

#endif

