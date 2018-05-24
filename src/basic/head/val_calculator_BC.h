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


#ifndef __VAL_CALCULATOR_BC_H
#define __VAL_CALCULATOR_BC_H

#include "val_calculator.h"

// #include "wall.h"

class Wall;

//---- newclass -----------------------------------------------------------

/*!
 * General parent-class for saving 
 * the pair-specific value of a quantity of the boundary particle used
 * for applying a boundary condition (BC) in each pair of particles.
 * Which kind of quantity and which kind of BC 
 * we are dealing with is specified by the child
 * classes.
 */
class ValCalculatorBC : public ValCalculatorPair
{
  protected:
    /*!
   * Are the BCs applied on free or frozen boundary particles?
     */
    bool m_frozen;

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
    vector<vector<pair<Wall*, double> >* > m_wallTable;

    /*!
     * Internal helper storing the \a Wall s in the periodic neighbouhood of each 
     * \a Boundary \a Particle and the distance to it
     */
    vector<vector<pair<Wall*, double> >* > m_periodicWallTable;
    
    /*!
     * Initialise the property list
     */
    virtual void init();
  
    /*!
     * Helper function for returning a copy of itself
     */
    virtual ValCalculatorPair* copyMySelf() = 0;
    
    /*!
     * Find the wall-segment where the rij-vector is hitting and the distance of ri to it
     * @param pD \a Pairdist for which to compute the distances
     * @param innerDist variable for distance of inner particle to wall
     * @param outerDist variable for distance of wall particle to wall
     */ 
    virtual void computeDists(Pairdist* pair, double& innerDist, double& outerDist); //FIXME: inline ?
    
    /*!
     * Adds the expressions used by this \a Symbol to the given list. 
     * This is only for symbols used in runtime-compiled expressions, so none here
     * @param usedSymbols List to be filled with own instances of \a TypedValue
     */
    virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
    {
      
    }

    /*!
     * In addition to calling the same method of its parant class (see 
     * documentatio there), this method registers the copy \a vc. And 
     * runs assertion for its own members variables.
     * @param[in] vc Copy of this \a ValCalculatorBC to be checked
     */
    virtual void makeCopySafe(ValCalculatorPair* vc) const;
    
  public:
  /*!
   * Constructor
   * @param symbol The symbol name
   */
    ValCalculatorBC(string symbol);

  /*!
     * Constructor for Node hierarchy
   */
    ValCalculatorBC(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ValCalculatorBC() {
    }
    
    /*!
     * Register this calculator and save the offset for the data stored
     * in a \a Pairdist in the argument \a slot
     */
    void setSlot(ColourPair* cp, size_t& slot, bool oneProp);
    
  /*!
     * Returns "BC"
   */
    virtual string myName() {
      return "BC";
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();

    /*!
     * For each wall-particle, find the walls in range
     */
    virtual void setupAfterParticleCreation();
};

#endif
