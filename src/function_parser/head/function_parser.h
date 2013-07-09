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



#ifndef __FUNCTION_PARSER_H
#define __FUNCTION_PARSER_H

#include <set>
#include <list>
#include <string>
#include <fstream>

using namespace std;

#include "function_node.h"

#include "fp_scalar.h"
#include "typed_value.h"
#include "binary_operators.h"

#define FUNCTION_PARSER_LOG(where, what)  FunctionParser::s_logFileStream << "{" << where << "} " << what << endl
#define FUNCTION_PARSER_LOG_ERROR(where, what)  FunctionParser::s_logFileStream << "{" << where << "} ERROR: " << what << endl

typedef list<TypedValue*> typed_value_list_t;

class UnknownSymbolCallback
{
 public:
  UnknownSymbolCallback() {}
  virtual ~UnknownSymbolCallback() {}

  virtual TypedValue *operator()(const string &e) = 0;
};


/*!
 * The function parser class takes a function definition as a string
 * and converts it into an internal tree structure. Thus, the function
 * can be either evaluated directly (slow) or compiled externally
 * by outputting and compiling c-code.
 *
 * The parser class is type aware.
 */
class FunctionParser
{
 protected:
  /*!
   * List of symbols (variables)
   */
  typed_value_list_t m_symbols;

  /*!
   * List of symbols that have been used in the expression
   */
  typed_value_list_t m_used_symbols;

  /*!
   * Function that is called when an unknown symbol is
   * encountered.
   */
  UnknownSymbolCallback *m_callback;
  
  /*!
   * The starting parse node
   */
  FunctionNode *m_main_node;

  /*!
   * The original expression that has been parsed
   */
  string m_expression;


  /*!
   * Return a constant/variable from either the symbol list
   * or by determining this is a number
   */
  FunctionNode *valueFromString(string expression);

  /*!
   * Parse this partial expression (recursive)
   */
  FunctionNode *parseThis(string expression);

 public:
  /*!
   * Constructor
   */
  FunctionParser();

  /*!
   * Destructor
   */
  virtual ~FunctionParser();

  /*!
   * Set the callback function
   */
  void setCallback(UnknownSymbolCallback *callback) {
    if (m_callback)
      delete m_callback;

    m_callback = callback;
  }

  /*!
   * Parse this expression
   */
  void parse(string expression);

  /*!
   * Return the value
   */
  Variant value() const {
    return m_main_node->value();
  }

  /*!
   * Return the value (only if scalar)
   */
  double valueAsScalar() const {
    return m_main_node->value().scalar();
  }

  /*!
   * Output C code
   */
  Variant toC() const {
    
/*    Variant test = m_main_node->toC();
    MSG_DEBUG("FunctionParser::toC()", "first C-expression: " << test.strings()[0]);*/
    return m_main_node->toC();
  }

  /*!
   * Output C code (only if scalar)
   */
  string toCAsScalar() const {
    return m_main_node->toC().scalarString();
  }

  /*!
   * Return the original expression
   */
  string expression() const {
    return m_expression;
  }

  /*!
   * Add a symbol to the list of symbols
   */
  void addSymbol(TypedValue *sym) {
    m_symbols.push_back(sym);
  }

  /*!
   * Enumerate the symbol that have been used in the expression
   */
  const typed_value_list_t &usedSymbols() const {
    return m_used_symbols;
  }

  /*!
  * Clear the symbols
  */
  void clearSymbols();


  /*!
   * Stream for the log file. 
   * It is static because we want only one file and one stream for all parsers
   */
  static ofstream s_logFileStream;
  /*!
   * Stream for the log file. 
   * Must be a pointer because there are modules using the copy-constructor 
   * of \a FunctionParser, and an ofstream is non-copyable!
   */
  /*   ofstream* m_logFileStream; */
  
  /*!
   * Counts the number of \a FunctionParser s; needed for \a s_logFileStream
   */
  static size_t s_numFunctionParsers;
  

};

#endif
