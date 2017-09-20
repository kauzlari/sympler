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


#ifndef __SYMBOL_H
#define __SYMBOL_H

#include "node.h"
// #include "simulation.h"

#include "smart_enum.h"

#include "data_format.h"

// #include "simulation.h"

// class Node;

/*!
 * Expression, computetd and saved in a symbol, to be used in expression 
 * of other modules.
 * This class is mainly used as a common interface for the parsing hierarchy.
 * The computation of the symbols is rather based on the sudivision into 
 * \a ValCalculator s and \a ParticleCache s, which inherit from Symbol.
*/

class Simulation;

class Symbol : public Node
{
  protected:
  /*!
   * The stage parameter tells when to compute the symbol relative to other symbols. This is 
   * important due to Symbols, which may depend on the previous 
   * computation of other Symbols
   */
    /*size_t*/int m_stage;
    
  /*!
   * The phase parameter tells when to compute the symbol in the timestep.
   */
    /*size_t*/int m_phase;
    
  /*!
    * The user definition for \a m_phase.
   */
    /*size_t*/int m_phaseUser;
    
    /*!
    * The name of the computed symbol to be used in other expressions
    */
    string m_symbolName;
    
    /*!
    * The datatype of this Symbol (usually DataFormat::DOUBLE, POINT, TENSOR)
    */
    DataFormat::datatype_t m_datatype;

#ifdef _OPENMP
    /*!
     * The number of doubles that each \a ValCalculator calculates
     */
//     int m_num_doubles;
#endif    

    /*!
   * Is this \a Symbol allowed to overwrite already existing symbols 
   * with name \a m_symbolName ?
   * Note: In this class we do not yet let the user control this member by 
   * an XML-attribute but only in selected sub-classes
     */
    bool m_overwrite;

    /*!
     * Initialise the PropertyList.
     */
    void init();

    /*!
     * Takes the \a name of one \a Symbol running at user stage 1 (\a m_phaseUser = 1 or 2) 
     * and tries to update \a m_stage of this \a ValCalculator
     */
    virtual void findStageForSymbolName(string name, bool& tooEarly, bool& nothing);

    /*!
     * Takes the \a name of one \a Symbol running at user stage 0 (\a m_phaseUser = 0 or 2) 
     * and tries to update \a m_stage of this \a ValCalculator
     */
    virtual void findStageForSymbolName_0(string name, bool& tooEarly, bool& nothing);

    
  public:
  /*!
   * Constructor for \a Node hierarchy
   * @param parent The parent node
   */
    Symbol(/*Node*/Simulation* parent);
   
  /*!
     * Constructor
   */
    Symbol(string symbol);
   
  /*!
     * Destructor
   */
    virtual ~Symbol() {
    }
  
  /*!
     * Setup this Symbol
   */
    virtual void setup() = 0;
    
    /*!
     * Return a string identifier for the calculator of this symbol
   */
    virtual string myName()
    {
      return m_properties.className(); 
    }

    /*!
     * Return the name of the computed symbol to be used in other expressions
     */
    virtual string mySymbolName()
    { 
      return m_symbolName;
    }
#ifdef _OPENMP
    /*!
     * Return the number of doubles that each \a ValCalculator calculates.
     * (Used by the ValCalculators themselves)
     */
//      virtual int numberOfDoubles() {
//        return m_num_doubles;
//      }

    /*!
     * Returns the datatype of this symbol.
     */
     virtual DataFormat::datatype_t symbolType() {
       return m_datatype;
     }

    /*!
     * Gets the number of doubles for a datatype, set by the \a DataFormat
     */
     virtual int setNumOfDoubles();
#endif

    /*!
     * Return the name of the computed symbols to be used in other expressions. This 
     * function takes into account that children of this class could compute more than one symbol
     */
    virtual list<string> mySymbolNames()
    {
      list<string> temp;
      assert(temp.empty());
      temp.push_back(m_symbolName);
      return temp;
    }

  /*!
     * Return the stage at which this calculator is being called.
   */
    virtual int stage() {
      return m_stage;
    }
    
    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    virtual bool findStage()
    {
      return true;
    }
    
    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    virtual bool findStage_0()
    {
      return true;
    }

    /*!
     * Checks for a consistent setup. Currently (2009-08-10) 
     * just needed for the tripletCalculator 
     */
    virtual void checkConsistency() {

    }

    virtual const bool doesOverwrite() const {
      return m_overwrite;
    }

};

/*--- Factory --- */

class SymbolFactory: public SmartEnum<SymbolFactory>
{
  public:
    virtual Symbol *instantiate(Simulation *parent) const = 0;

  protected:
    SymbolFactory(const string &name)
  : SmartEnum<SymbolFactory>(name) { }
};


template <class T>
    class SymbolRegister: public SymbolFactory
{
  public:
    SymbolRegister(const string &name)
  : SymbolFactory(name) { }

    virtual Symbol *instantiate(Simulation *parent) const {
      return new T(parent);
    }
};

#endif
