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



#ifndef __FUNCTION_COMPILER_H
#define __FUNCTION_COMPILER_H 

#include <set>
#include <vector>
#include <ostream>
#include <functional>

#include <dlfcn.h>
#include <stdlib.h>

#include "function_parser.h"

using namespace std;

#define C_FC_FN_NAME "fn"
#define C_FC_PREFIX "__function_compiler_tmp_"

typedef void (*gfnptr_t)(...);

/*!
 * Externally compiles a function using an appropriate C compiler
 * and loads the created object file.
 */
class FunctionCompiler
{
 protected:
  /*!
   * How many function objects are there (to keep track of .c files)
   */
  static int s_function_counter;   

  /*!
   * A pointer to the parser object holding the actual function
   */
  FunctionParser *m_parser;

  /*!
   * The header of the external function
   */
  string m_header;

  /*!
   * The variables the result should be stored in
   */
  vector<string> m_result_vars;

  /*!
   * Name of the shared object 
   * (to be deleted when Function is destroyed).
   */
  string m_so_filename;

  /*!
   * Handle to the compiled shared object.
   */
  void *m_handle;

  /*!
   * Pointer to the function inside the shared object.
   */
  gfnptr_t m_fn;

  /*!
   * Compile the function externally
   */
  bool compile();

  /*!
   * Closes the file handle.
   */
  void close();

public:
  /*! 
   * Constructor
   */
  FunctionCompiler();

  /*!
   * Destructor
   */
  ~FunctionCompiler();

  /*!
   * Use this function if only one result (i.e. no vector) is
   * expected.
   * Note: This is not checked for consistency!!!
   */
  void setResultString(string result_var) {
    m_result_vars.clear();
    m_result_vars.push_back(result_var);
  }

  /*!
   * Set the result string for the external function.
   * Note: This is not checked for consistency!!!
   */
  void setResultStrings(vector<string> result_vars) {
    m_result_vars = result_vars;
  }

  /*!
   * Set the header string for the external function.
   * Note: This is not checked for consistency!!!
   */
  void setHeader(string header) {
    m_header = header;
  }

  /*!
   * Set the actual expression and compile
   */
  void setParserAndCompile(FunctionParser *parser);

  gfnptr_t fn() const {
    return m_fn;
  }

  bool isNull() const {
    return !m_handle;
  }
};


/*
template<typename T>
ostream& operator<<(ostream &out, const FunctionCompiler f) {
  out << f.expression();
	
  return out;
}
*/

#endif
