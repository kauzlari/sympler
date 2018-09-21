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


#ifndef __SYMBOL_H
#define __SYMBOL_H

#include "node.h"
#include "data_format.h"
#include "function_parser.h"


// class Node;

/*!
 * Expression, computed and saved in a symbol, to be used in expressions 
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
  int m_stage;
  
  /*!
   * The phase parameter tells when to compute the symbol in the timestep.
   */
  int m_phase;
  
  /*!
   * The user definition for \a m_phase.
   */
  int m_phaseUser;
  
  /*!
   * The name of the computed symbol to be used in other expressions
   */
  string m_symbolName;
  
  /*!
   * The datatype of this Symbol (usually DataFormat::DOUBLE, POINT, TENSOR)
   */
  DataFormat::datatype_t m_datatype;
  
  /*!
   * Is this \a Symbol allowed to overwrite already existing symbols 
   * with name \a m_symbolName ?
   * Note: In this class we do not yet let the user control this member by 
   * an XML-attribute but only in selected sub-classes
   */
  bool m_overwrite;

  /*!
   * Will the symbol introduced by this \a Symbol be protected from 
   * automatic reset (=0)? The default value in the constructor is 
   * 'false' for class \a Symbol.
   * FIXME: currently (2018-05-09) not used, but could be useful, 
   * especially if we add the feature, that the \a Symbol takes itself 
   * control over the resetting while \a m_persistency = true to 
   * protect from automatic reset. First incomplete code related to 
   * this feature exists in some \a Symbol s, but is commented out due 
   * to a vanished need for it. 
   */
  bool m_persistency;
  
  /*!
   * This string holds the symbols, which are not waited for to be computed beforehand
   */
  string m_oldSymbols;
  
  /*!
   * Initialise the PropertyList.
   */
  void init();
  
  /*!
   * Deals with setting of \a m_stage at user stage 1 (\a m_phaseUser = 1 or 2)
   * depending on value of \a m_overwrite
   * @param tooEarly Modified to 'true' if stage cannot be determined yet
   * @param nothing Set to false if other \a Symbol found, which computes \a mySymbolNames()
   */
  virtual void checkOverwriteForStageFinding(bool& tooEarly, bool& nothing);
  
  /*!
   * Deals with setting of \a m_stage at user stage 0 (\a m_phaseUser = 0 or 2)
   * depending on value of \a m_overwrite
   * FIXME: Try to avoid code duplication for _0 versions
   * @param tooEarly Modified to 'true' if stage cannot be determined yet
   * @param nothing Set to false if other \a Symbol found, which computes \a mySymbolNames()
   */
  virtual void checkOverwriteForStageFinding_0(bool& tooEarly, bool& nothing);
  
  /*!
   * Takes the \a name of one \a Symbol running at user stage 1 (\a m_phaseUser = 1 or 2) 
   * and tries to update \a m_stage of this \a ValCalculator
   * @param name Name of computed symbol to be checked for stage of computation
   * @param tooEarly Modified to 'true' if stage cannot be determined yet
   * @param nothing Set to false if other \a Symbol found, which computes \a name
   */
  virtual void findStageForSymbolName(string name, bool& tooEarly, bool& nothing);
  
  /*!
   * Takes the \a name of one \a Symbol running at user stage 0 (\a m_phaseUser = 0 or 2) 
   * and tries to update \a m_stage of this \a ValCalculator
   * FIXME: Try to avoid code duplication for _0 versions
   * @param name Name of computed symbol to be checked for stage of computation
   * @param tooEarly Modified to 'true' if stage cannot be determined yet
   * @param nothing Set to false if other \a Symbol found, which computes \a name
   */
  virtual void findStageForSymbolName_0(string name, bool& tooEarly, bool& nothing);
  
  /*!
   * Adds the expressions used by this \a Symbol to the given list. 
   * There are subclasses, which indeed do not leave the list in its original state
   * @param usedSymbols List to be filled with own instances of \a TypedValue
   * FIXME: this construction is not optimal, since only those subclasses which 
   * use it should have to know about the type typed_value_list_t
   * FIXME: When all inheriting classes are transferred to Symbol::findStage/-_0(), 
   * then remove exception and make this method pure virtual, such that the programmer 
   * must think about the appropriate behaviour for every child 
   */
  virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
  {
    throw gError("Symbol::addMyUsedSymbolsTo", "Error for module " + className() + ": This method should not have been called. This means that a programmer forgot to implement this method for an inheriting class. In the simplest use case the missing implementation would just do nothing. If you have no idea what to do, find the person to blame for the affected class/module. Oh, it's you, well...");
  }

  /*!
   * Returns the strings of those \a Symbols that the given class depends on
   * due to hard-coded reasons (not due to runtime compiled expressions).
   */
  virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
  {
    throw gError("Symbol::myHardCodedDependencies", "Error for module " + className() + ": This method should not have been called. This means that a programmer forgot to implement this method for an inheriting class. In the simplest use case the missing implementation would just do nothing. If you have no idea what to do, find the person to blame for the affected class/module. Oh, it's you, well...");
  }

  /*!
   * The returned string contains those terms from runtime compiled expressions, 
   * which should be ignored when determining the stage. The expressions are separated by " | "
   * An "empty" string must have the form "---".
   */
  virtual string usedSymbolsIgnoredForStaging() const {

    if(m_oldSymbols == "---")
      return m_symbolName;
    else return string(m_oldSymbols + "|" + m_symbolName);
  }
  
  /*!
   * Helper function which removes indices and brackets from single terms in \a Function 
   * expressions. Subclasses must decide which kind of cleaning to perform
   * @param name Single term from a \a Function expression
   */
  virtual void cleanSymbol(string& name) const {
  }

  // see CONVENTION5 for rule about persistencies
  /*
   * If this is true, the symbol cannot be cleared
   * Currently (2018-02-20) = false for all \a Symbol s 
   * (see symbol.cpp)
   */
  static const bool s_persistency;

  
 public:
  
  /*!
   * Constructor for \a Node hierarchy
   * @param parent The parent node
   */
  Symbol(Simulation* parent);
  
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
   * Initialise all variables of this \a Symbol
   * FIXME: This function was added pretty recently (2018-02-19), and 
   * some of the children do not call it in their own setup(). It also
   * does not yet setup all members it declares.
   */
  virtual void setup();
  
  /*!
   * Return a string identifier for the calculator of this symbol
   */
  virtual string myName() const
  {
    return m_properties.className(); 
  }
  
  /*!
   * Return the name of the computed symbol to be used in other expressions
   */
  virtual string mySymbolName() const
  { 
    return m_symbolName;
  }
  
  /*!
   * Lets clients check if an instance is overwriting
   */
  virtual const bool& overwriting() const {
    return m_overwrite;
  }

  /*!
   * Returns the datatype of this symbol.
   */
  virtual DataFormat::datatype_t symbolType() {
    return m_datatype;
  }
  
#ifdef _OPENMP

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
   * Implements the main stage finding logic for user-stage 1. Calls methods
   * that may be overridden by child classes.
   */
  virtual bool findStage();
  
  /*!
   * Determines \a m_stage of the current \a Symbol.
   * Implements the main stage finding logic for user-stage 0. Calls methods
   * that may be overridden by child classes.
   * FIXME: Try to avoid code duplication for _0 versions
   */
  virtual bool findStage_0();
  
  /*!
   * Checks for a consistent setup. Currently (2009-08-10) 
   * just needed for the \a TripletCalculator 
   */
  virtual void checkConsistency() {
    
  }
  
  virtual const bool doesOverwrite() const {
    return m_overwrite;
  }

  /*!
   * Helper to remove string names of symbols given by 
   * \a symbolChainToBeRemoved  from a list of symbols 
   * \a symbolsFromWhichToRemove.
   * @param[in] symbolChainToBeRemoved The symbols to be removed, given 
   * as a string with separator "|" between the names of the symbols
   * @param[out] symbolsFromWhichToRemove Contains the list of symbols 
   * from which some of them should be removed
   */
  static void removeFromSymbolList(string symbolChainToBeRemoved, list<string>& symbolsFromWhichToRemove);
  
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
