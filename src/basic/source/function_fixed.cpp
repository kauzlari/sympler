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


#include "function_fixed.h"

#include "fp_scalar.h"

#include <algorithm>

extern bool fexists(const char *filename);



/* FunctionFixed */

FunctionFixed::FunctionFixed()
{
}


FunctionFixed::~FunctionFixed()
{
}


/* Methods */

int FunctionFixed::addVariable(string name)
{
//   MSG_DEBUG("FunctionFixed::addVariable", "m_expression = " << m_expression);
  m_vars.push_back(name);
//   MSG_DEBUG("FunctionFixed::addVariable", "added variable " << name << ", size now = " << m_vars.size());

  return m_vars.size()-1;
}

int FunctionFixed::addVariables(string name1, string name2)
{
  m_vars.push_back(name1);
  m_vars.push_back(name2);

  return m_vars.size()-1;
}

int FunctionFixed::addVariables(string name1, string name2, string name3)
{
  m_vars.push_back(name1);
  m_vars.push_back(name2);
  m_vars.push_back(name3);

  return m_vars.size()-1;
}


void FunctionFixed::deleteVariable(string name)
{
  vector<string>::iterator i;

  i = find(m_vars.begin(), m_vars.end(), name);
  if (i != m_vars.end())
    m_vars.erase(i);
  else
    throw gError("FuntionFixed::deleteVariable", "Variable '" + name + "' not found.");
}


void FunctionFixed::compile()
{
  string header = "";

//   m_expression = expression;

  m_used_vars.clear();

  if (m_vars.empty())
    m_vars.push_back("x");

  header = "double *result";
  
  // if this is a recompilation, then make sure, the symbol lists are clean
//   m_parser.clearSymbols();
  
  for (vector<string>::iterator i = m_vars.begin(); i != m_vars.end(); ++i) {
    header += ", double " + *i;

    m_parser.addSymbol(new FPScalarVariable(*i, NULL));
  }
  m_parser.parse(m_expression);

  for (typed_value_list_t::const_iterator i = m_parser.usedSymbols().begin();
       i != m_parser.usedSymbols().end(); ++i)
    m_used_vars.insert((*i)->name());

  m_compiler.setResultString("*result");
  m_compiler.setHeader(header);
  m_compiler.setParserAndCompile(&m_parser);
}

