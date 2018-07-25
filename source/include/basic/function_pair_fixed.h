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



#ifndef __FUNCTION_PAIR_FIXED_H
#define __FUNCTION_PAIR_FIXED_H 

#include "function_pair.h"

/*!
 * A function which takes a pair as its argument and, additionally, one argument with a meaning predefined by the client object
 */
class FunctionPairFixed: public FunctionPair
{
  protected:
  /*!
   * Fixed variables that can be used.
   */
    vector<string> m_vars;           

  /*!
     * Fixed variables that are actually being used.
   */
    set<string> m_used_vars;


  public:
  /*!
   * Constructor
   */
    FunctionPairFixed();

  /*!
     * Destructor
   */
    virtual ~FunctionPairFixed();

  /*!
     * Compile the function
   */
    virtual void compile();

    /*!
     * Sets only the expression, without adding this \a Function to the 
     * \a Function::toCompile list
     * @param expr The expression
     */
    void setOnlyExpr(string expr)
    {
      m_expression = expr; 
    }
  
    /*!
     * Add a variable of name \a name
     * @param name Name of the variable
     */
    int addVariable(string name);

    /*!
     * Delete variable with name \a name
     * @param name Name of the variable
     */
    void deleteVariable(string name);

    /*!
     * Return a set of all variables actually used in the expression
     */
    const set<string> &usedVariables() const {
      return m_used_vars;
    }
    
    /*!
     * Call the function
     * @param val Return value
     * @param pr Pairdist this function takes as the argument
   */
    void operator()(void* val, double* x, Pairdist* pr) const {
      m_compiler.fn()
          (val,
           x,
           &(pr->m_distance),
           pr->tag.data(),
           pr->firstPart(),
           pr->secondPart(),
           pr->firstPart()->tag.data(),
           pr->secondPart()->tag.data());
    }
};


template<typename T>
    ostream& operator<<(ostream &out, const FunctionPairFixed &f) {
  out << f.expression();
	
  return out;
    }



#endif
