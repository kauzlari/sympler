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



#ifndef __FUNCTION_H
#define __FUNCTION_H

#include <string>

#include "function_parser.h"
#include "function_compiler.h"

using namespace std;

/*!
 * Base class for user-defined functions
 */
class Function
{
protected:
  /*!
   * The expression to compile
   */
  string m_expression;

  /*!
   * The function parser
   */
  FunctionParser m_parser;

  /*!
   * The function compiler
   */
  FunctionCompiler m_compiler;
                                        
public:
  /*!
   * Constructor
   */
  Function();

  /*!
   * Destructor
   */
  virtual ~Function();
  
  /*!
   * Compile this function using gcc
   */
  virtual void compile() = 0;

  /*!
   * Set the expression
   */
  void setExpression(string expression);

  /*!
   * Return the (unparsed) expression
   */
  const string &expression() const {
    return m_expression;
  }

  /*!
   * Is the expression null or not set?
   */
  bool isNull() const {
    return m_expression == "";
  }

    /*!
   * Enumerate the symbols that have been used in this \a Function
     */
  const typed_value_list_t &usedSymbols() const {
    return m_parser.usedSymbols();
  }

  
  /*!
   * List of functions to be compiled by the \a Controller
   */
  static set<Function*> toCompile;
};

#endif
