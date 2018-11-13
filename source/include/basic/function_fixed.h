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



#ifndef __FUNCTION_FIXED_H
#define __FUNCTION_FIXED_H 

#include "function.h"
#include "function_parser.h"
#include "function_compiler.h"


/*!
 * A function with a fixed number of double arguments
 * that returns a double.
 */
class FunctionFixed : public Function
{
 protected:
  /*!
   * Variables that can be used.
   */
  vector<string> m_vars;           

  /*!
   * Variables that are actually being used.
   */
  set<string> m_used_vars;

 public:
  /*!
   * Constructor
   */
  FunctionFixed();

  /*!
   * Destructor
   */
  ~FunctionFixed();

  /*!
   * Compile this function to C-code
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
   * Add two variables of name \a name1 and \a name2
   * @param name1 Name of the first variable
   * @param name2 Name of the second variable
   */
  int addVariables(string name1, string name2);

  /*!
   * Add three variables of names \a name1, \a name2 and \a name3
   * @param name1 Name of the first variable
   * @param name2 Name of the second variable
   * @param name3 Name of the third variable
   */
  int addVariables(string name1, string name2, string name3);

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
   * Call the compiled function (one argument)
   * @param x First argument
   */
  double operator()(double x) const {
    assert(m_vars.size() == 1);

    double result;
    m_compiler.fn()(&result, x);
    return result;
  }

  /*!
   * Call the compiled function (two arguments)
   * @param x First argument
   * @param y Second argument
   */
  double operator()(double x, double y) const {
    assert(m_vars.size() == 2);

    double result;
    m_compiler.fn()(&result, x, y);
    return result;
  }

  /*!
   * Call the compiled function (three arguments)
   * @param x First argument
   * @param y Second argument
   * @param z Third argument
   */
  double operator()(double x, double y, double z) const {
    assert(m_vars.size() == 3); 

    double result;
    m_compiler.fn()(&result, x, y, z);
    return result;
  }

  /*!
   * Call the compiled function (four arguments)
   * @param x First argument
   * @param y Second argument
   * @param z Third argument
   * @param t Fourth argument
   */
  double operator()(double x, double y, double z, double t) const {
    assert(m_vars.size() == 4);

    double result;
    m_compiler.fn()(&result, x, y, z, t);
    return result;
  }

  /*!
   * Call the compiled function (five arguments)
   * @param x First argument
   * @param y Second argument
   * @param z Third argument
   * @param t Fourth argument
   * @param u Fifth argument
   */
  double operator()(double x, double y, double z, double t, double u) const {
    assert(m_vars.size() == 5);

    double result;
    m_compiler.fn()(&result, x, y, z, t, u);
    return result;
  }

  /*!
   * Call the compiled function (six arguments)
   * @param x First argument
   * @param y Second argument
   * @param z Third argument
   * @param t Fourth argument
   * @param u Fifth argument
   * @param v Sixth argument
   */
  double operator()(double x, double y, double z, double t, double u, double v) const {
    assert(m_vars.size() == 6);

    double result;
    m_compiler.fn()(&result, x, y, z, t, u, v);
    return result;
  }

  /*!
   * Call the compiled function (seven arguments)
   * @param x First argument
   * @param y Second argument
   * @param z Third argument
   * @param u Fourth argument
   * @param v Fifth argument
   * @param w Sixth argument
   * @param t Seventh argument
   */
  double operator()(double x, double y, double z, double u, double v, double w, double t) const {
    assert(m_vars.size() == 7);

    double result;
    m_compiler.fn()(&result, x, y, z, u, v, w, t);
    return result;
  }

//  /*!
//   * Has the function been compiled?
  //   */
//  bool isNull() const {
//    return m_compiler.isNull();
//  }
};


template<typename T>
ostream& operator<<(ostream &out, const FunctionFixed &f) {
  out << f.expression();
	
  return out;
}

#endif
