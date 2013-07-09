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


#include <stdio.h>

#include "function_parser.h"
#include "function_compiler.h"

#include "fp_scalar.h"
#include "fp_vector.h"

typedef double (*myfn_t)(double, double, double);

int main(int argc, char *argv)
{
  FunctionParser p;

  try {
    point_t x = {{{ 0, 0, 0 }}};
    string e = "x:sin(x)";

    p.addSymbol(new FPVectorVariable("x", &x));

    p.parse(e);

    cout << e << endl;
    cout << p.toCAsScalar() << endl;

    cout << "Interpreting" << endl;

    x.x = 0;
    cout << "x = " << x << " => f(x) = " << p.valueAsScalar() << endl;
    x.x = 2; x.y = 1;
    cout << "x = " << x << " => f(x) = " << p.valueAsScalar() << endl;
    x.x = 4;
    cout << "x = " << x << " => f(x) = " << p.valueAsScalar() << endl;

    cout << "Compiling" << endl;

    FunctionCompiler c;
    
    c.setHeader("double x0, double x1, double x2");
    c.setParserAndCompile(&p);

    x.x = 0; x.y = 0;
    cout << "x = " << x << " => f(x) = " << ((myfn_t) c.fn())(x.x, x.y, x.z) << endl;
    x.x = 2; x.y = 1;
    cout << "x = " << x << " => f(x) = " << ((myfn_t) c.fn())(x.x, x.y, x.z) << endl;
    x.x = 4;
    cout << "x = " << x << " => f(x) = " << ((myfn_t) c.fn())(x.x, x.y, x.z) << endl;
  } catch (gError &err) {
    cout << "The following error occured: " << endl << err.message() << endl;
  }
}
