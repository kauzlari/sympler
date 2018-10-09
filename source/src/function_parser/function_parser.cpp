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


#include <algorithm>

using namespace std;

#include "function_parser.h"

#include "fp_scalar.h"
#include "unary_functions.h"
#include "unary_operators.h"
#include "binary_operators.h"

size_t FunctionParser::s_numFunctionParsers = 0;
ofstream FunctionParser::s_logFileStream;

FunctionParser::FunctionParser(): m_callback(NULL), m_main_node(NULL)
{
  ++FunctionParser::s_numFunctionParsers;

  // If no other FunctionParser has opened the static log file stream, do it
  if(!FunctionParser::s_logFileStream.is_open()) {
    // assert for consistency-check
    assert(FunctionParser::s_numFunctionParsers == 1);
    FunctionParser::s_logFileStream.open("function_parser.log");
    FunctionParser::s_logFileStream.flags(ios::scientific);
  }
}


FunctionParser::~FunctionParser()
{
  --FunctionParser::s_numFunctionParsers;

  // Am I the last FunctionParser to be destroyed?
  if(FunctionParser::s_numFunctionParsers == 0) {
    // consistency-check since with the current (2013-06-07) implementation
    // the stream should still be open
    assert(FunctionParser::s_logFileStream.is_open());
    // then close log file stream if still open
    FunctionParser::s_logFileStream.close();    
  }
 
  if (m_callback)
    delete m_callback;

  if (m_main_node)
  {
    // symbols may NOT be deleted here, because still accessed in m_symbols
    if(!dynamic_cast<TypedValue*>(m_main_node)) delete m_main_node;
  }
  
  for (list<TypedValue*>::iterator i = m_symbols.begin(); i != m_symbols.end(); ++i)
  {
    delete *i;
  }

}

void FunctionParser::clearSymbols()
{
  if (m_main_node)
    // symbols may NOT be deleted here, because still accessed in m_symbols
    if(!dynamic_cast<TypedValue*>(m_main_node)) delete m_main_node;
  
  
  for (list<TypedValue*>::iterator i = m_symbols.begin(); i != m_symbols.end(); ++i)
  {
    delete *i;
  }
  m_main_node = NULL;
  m_symbols.clear();
  // deleting in m_used_symbols is not necessary because there are only 
  // pointers to entries in m_symbols
  m_used_symbols.clear();
}

void FunctionParser::parse(string expression)
{
//   FUNCTION_PARSER_LOG("FunctionParser::parse", "START: expression = " << expression);
  if (m_main_node)
  {
/*    FUNCTION_PARSER_LOG("FunctionParser::parse", "DELETING main node, BEFORE m_symbols = ");
    for (list<TypedValue*>::iterator i = m_symbols.begin(); i != m_symbols.end(); ++i)
    {
      FUNCTION_PARSER_LOG("FunctionParser::parse", (**i).m_name);
    }*/
    
    // symbols may NOT be deleted here, because still accessed in m_symbols
    if(!dynamic_cast<TypedValue*>(m_main_node)) delete m_main_node;
    
/*    FUNCTION_PARSER_LOG("FunctionParser::parse", "DELETED main node, NOW m_symbols = ");
    for (list<TypedValue*>::iterator i = m_symbols.begin(); i != m_symbols.end(); ++i)
    {
      FUNCTION_PARSER_LOG("FunctionParser::parse", (**i).m_name);
    }*/
  
  }
    
  m_used_symbols.clear();
  m_expression = expression;

  try {
//     FUNCTION_PARSER_LOG("FunctionParser::parse", "PARSING main node");
    m_main_node = parseThis(expression);
  } catch (gError &err) {
    FUNCTION_PARSER_LOG_ERROR
      ("FunctionParser::parse",
       "Error while parsing '" + expression + "'. " +
       "Subroutine: " + err.where() + ", Message: " + err.msg());
    throw gError
      ("FunctionParser::parse",
       "Error while parsing '" + expression + "'. " +
       "Subroutine: " + err.where() + ", Message: " + err.msg());
  }
//   FUNCTION_PARSER_LOG("FunctionParser::parse", "END: expression = " << expression);
}


/*size_t*/int findWithoutParentheses(const string s, const string substr)
{
//    FUNCTION_PARSER_LOG("FunctionParser::findWithoutParentheses", "START: s = " << s << ", substr = " << substr);
  // new: we search from right to left
  int p = s.size()-1;
  size_t level = 0;

  while ((p > -1) && (string(s, p, substr.size()) != substr || level != 0) /*&& p > -1*/) {
    if (s[p] == ')')
      level++;
    else if (s[p] == '(')
      level--;

    --p;
  }
  
  // old
/*  size_t p = 0;
  size_t level = 0;

  while ((string(s, p, substr.size()) != substr || level != 0) && p < s.size()) {
    if (s[p] == '(')
      level++;
    else if (s[p] == ')')
      level--;

    p++;
  }*/

//      FUNCTION_PARSER_LOG("FunctionParser::findWithoutParentheses", "found p = " << p);
  return p;
}


struct find_symbol: public unary_function<const TypedValue *, bool> {
  string m_expr;
  find_symbol(const string &expr): m_expr(expr) { }
  bool operator()(const TypedValue *v) {
    return v->name() == m_expr;
  }
};

FunctionNode *FunctionParser::valueFromString(string expression)
{
//    FUNCTION_PARSER_LOG("FunctionParser::valueFromString", "START: expression = " << expression << endl << "My symbols:");
//    for (list<TypedValue*>::iterator i = m_symbols.begin(); i != m_symbols.end(); ++i)
//   {
//     cout << (*i)->name() << endl;
//   }

  
  list<TypedValue*>::iterator sym =
    find_if(m_symbols.begin(), m_symbols.end(), find_symbol(expression));
  if (sym == m_symbols.end()) {
    /* Okay, this is not in our list of symbols, try to request one */

    TypedValue *n = NULL;
    char *endptr = NULL;
    double d;

    if (m_callback) {
      n = (*m_callback)(expression);

      if (n) {
//         FUNCTION_PARSER_LOG("FunctionParser::valueFromString", "callback-case: pushing and returning symbol " << n->m_name);
  /* So that we'll find this next time. */
	m_symbols.push_back(n);
	
	return n;
      }
    }
    
    // try to convert the string to a number
    // strtod would convert "" to "0"
    if(expression == "") {
      FUNCTION_PARSER_LOG_ERROR("FunctionParser::valueFromString", "I have to parse the empty expression \"\". Did you forget an argument of a function or operator?");
      throw gError("FunctionParser::valueFromString", "I have to parse the empty expression \"\". Did you forget an argument of a function or operator?");
    }    

//     FUNCTION_PARSER_LOG("FunctionParser::valueFromString", "Converting " << expression.c_str() << " to a double");   
    
    d = strtod(expression.c_str(), &endptr);
    
    string takeExpr("(" + expression + ")");

    if ((size_t) (endptr-expression.c_str()) != expression.size()) {
    // check for pi = 3.1415...
      if(expression == "pi")
      {
        d = M_PI;
        // this ensures the maximum number of digits
        takeExpr = "(M_PI)";
      }
      else
      {
        string symbol_list("");

        for (typed_value_list_t::iterator i = m_symbols.begin();
	        i != m_symbols.end(); ++i) {
	         if (i != m_symbols.begin()) symbol_list += ", ";
         symbol_list += "'" + (*i)->name() + "'";
        }
	// this will be catched and then this message also written with FUNCTION_PARSER_LOG
        throw gError
          ("FunctionParser::valueFromString",
          "'" + expression + "' is neither a defined symbol nor a number. "
          "Check whether the respective species support this expression. "
          "Unbalanced brackets might also be the reason. "
          "Currently initialised symbols are: " + symbol_list);
      }
    }
     FUNCTION_PARSER_LOG("FunctionParser::valueFromString", "constant-case. d = " << d << ", takeExpr = " << takeExpr << " for expression = " << expression);

     return new FPScalarConstant(takeExpr, d);
  } else {
    /* Okay, it's a symbol */

    if (find(m_used_symbols.begin(), m_used_symbols.end(), *sym) == m_used_symbols.end())
      m_used_symbols.push_back(*sym);
     FUNCTION_PARSER_LOG("FunctionParser::valueFromString", "symbol-case: returning symbol " << (*sym)->m_name);

    return *sym;
  }
}

/*bla*/
FunctionNode *FunctionParser::parseThis(string expression)
{
FUNCTION_PARSER_LOG("FunctionParser::parseThis", "start-expression = " << expression);
  string expr = expression;
  FunctionNode_Factory *factory = NULL;
  /*size_t*/int pos;

//   if (expr[0] == '(' && expr[expr.size()-1] == ')' && 
//    expr.find(")") == expr.size()-1)
//     expr = string(expr, 1, expr.size()-2);

  // is there an all enclosing bracket? we use 'do while', because some user could
  // write unnecessarily many brackets
  if(expr[0] == '(')
  {
    if(expr[1] == ')') {
      // this will be catched and then this message also written with FUNCTION_PARSER_LOG
      throw gError("FunctionParser::parseThis", "Empty bracket!");
    }
    bool foundBrackets;
    do
    {
      FUNCTION_PARSER_LOG("FunctionParser::parseThis", "in bracket loop with expression = " << expr);
      foundBrackets = false;
      size_t open = 1;
      size_t position = 0;
      // find the next bracket, either "(" or ")"
      while(open)
      {
        size_t another = 1+position;
        position = expr.find(")", another);
// FUNCTION_PARSER_LOG("FunctionParser::parseThis", "position = " << position);
        another = expr.find("(", another);
// FUNCTION_PARSER_LOG("FunctionParser::parseThis", "another = " << another);
        if(position > another)
        {
          position = another;
          ++open;
        }
        else --open;
// FUNCTION_PARSER_LOG("FunctionParser::parseThis", "open = " << open);
      }
      // now, the first bracket should be closed; if we are at the end, we have an 
      // all enclosing bracket that has to be removed
      if(position == expr.size()-1)
      {
        expr = string(expr, 1, expr.size()-2);
        foundBrackets = true;
      }
    } while (foundBrackets);
  }
  FUNCTION_PARSER_LOG("FunctionParser::parseThis", "after brackets check = " << expr);

  for (int i = 0; i < C_MAX_PRIORITY && !factory; ++i) {
    list<FunctionNode_Factory*> cur = FunctionNode_Factory::byPriority(i);

    for (list<FunctionNode_Factory*>::iterator cf = cur.begin();
	   cf != cur.end() && !factory; ++cf) {

      pos = findWithoutParentheses(expr, (*cf)->name());

      // new: we search from right to left
      if (pos != -1) 
        // old
//       if (pos != expr.size()) 
      {
        if((*cf)->nodeType() == FunctionNode::UNARY_FUNCTION)
        {
          // only if the unary function was found at the beginning of the expression, it is correct
          if(pos == 0) 
          {
            factory = *cf; 
//             FUNCTION_PARSER_LOG("FunctionParser::parseThis", "found factory " << (*cf)->name() << " with pos = " << pos);
          }
        }
        else
        {
//           FUNCTION_PARSER_LOG("FunctionParser::parseThis", "found factory " << (*cf)->name() << " with pos = " << pos);
	        factory = *cf;
        }
      }
    }
  }

  if (!factory) {
    /* No factory? This must be a variable or a constant. */
    FUNCTION_PARSER_LOG("FunctionParser::parseThis", "NO-FACTORY case");
    FunctionNode *node = valueFromString(expr);

    return node;
  } else {
    FunctionNode *node = factory->instantiate();

    switch (factory->nodeType()) {
    case FunctionNode::UNARY_FUNCTION: {
      FNUnaryFunction *n = dynamic_cast<FNUnaryFunction*>(node);
      FUNCTION_PARSER_LOG("FunctionParser::parseThis", "UNARY_FUNCTION case");

      if (pos != 0)
      {
        FUNCTION_PARSER_LOG_ERROR("FunctionParser::parseThis", "Aborting because pos!=0; expr = " << expr);
	// will be catched and then written also with FUNCTION_PARSER_LOG_ERROR
	throw gError
	  ("FunctionParser::parseThis",
	   "Syntax error.");
      }
      if (!n) {
	// will be catched and then written also with FUNCTION_PARSER_LOG_ERROR
	throw gError
	  ("FunctionParser::parseThis",
	   "(Internal error) Unary function announced for symbol '" + factory->name() +
	   "' but retrieved object is derived from a different class.");
      }
//       FUNCTION_PARSER_LOG("FunctionParser::parseThis", "UNARY-CASE");
      n->setUnary
	(parseThis(string(expr, factory->name().size(), expr.size())));
    }
      break;
    case FunctionNode::BINARY_OPERATOR: {
      FNBinaryOperator *n = dynamic_cast<FNBinaryOperator*>(node);
      FUNCTION_PARSER_LOG("FunctionParser::parseThis", "BINARY-CASE; expr = " << expr);
      
      if (!n) {
	// will be catched and then written also with FUNCTION_PARSER_LOG_ERROR
	throw gError
	  ("FunctionParser::parseThis",
	   "(Internal error) Binary operator announced for symbol '" + factory->name() +
	   "' but retrieved object is derived from a different class.");
      }
      
      /* We have to check for the '-' operator. Any way to do this more cleverly,
	 i.e., to not explicitly check for '-' here? */
      
      if (factory->name() == "-" && pos == 0) {
	//	n->setBinary(parseThis("0"), parseThis(string(expr, 1, expr.size())));
	delete node;
	
	FNNegation *new_n = new FNNegation();
	new_n->setUnary(parseThis(string(expr, 1, expr.size())));
	node = new_n;
      } else {
	n->setBinary
	  (parseThis(string(expr, 0, pos)),
	   parseThis(string(expr, pos+factory->name().size(), expr.size())));
      }
    }
      break;
    default: {
      // will be catched and then also written with FUNCTION_PARSER_LOG_ERROR
      throw gError
	("FunctionParser::parseThis",
	 "(Internal error) Don't know how to handle node type " +
	 ObjToString(factory->nodeType()) + " for symbol " + factory->name() + ".");
    }
    }
  
    return node;
  }

  return NULL;
}


/*static*/ void FunctionParser::addToTypedValueList(const typed_value_list_t& newSymbols, typed_value_list_t& usedSymbols) {
      
  for(typed_value_list_t::const_iterator s = newSymbols.begin(); s != newSymbols.end(); ++s)
    usedSymbols.push_back(*s);

}


/*static*/ void FunctionParser::removeFromTypedValues(string toRemove, typed_value_list_t& usedSymbols) {

  FUNCTION_PARSER_LOG("FunctionParser::removeFromTypedValues", "START: toRemove = '" << toRemove << "'");
      
  // go through the symbols in 'toRemove' and remove those from usedSymbols 
  if (toRemove != "---") {
    bool run = true;
    string working = toRemove;
    while(run) {
      string cur;
      size_t pos = working.find('|');
      
      if (pos == string::npos) {
	run = false;
	cur = working;
      }
      else {
	cur = string(working, 0, pos);
	working = string(working, pos+1);
      }
      
      typed_value_list_t symbolsToRemove;
      // determine what to remove
      for(typed_value_list_t::const_iterator s = usedSymbols.begin(); s != usedSymbols.end(); ++s) {
	FUNCTION_PARSER_LOG("FunctionParser::removeFromTypedValues", "now comparing used symbol '" << (*s)->name() << "' with to be removed symbol '" << cur << "'");
	if((*s)->name() == cur) {
	  symbolsToRemove.push_back(*s);
	  FUNCTION_PARSER_LOG("FunctionParser::removeFromTypedValues", "will remove" << (*s)->name());
	}
	else
	  FUNCTION_PARSER_LOG("FunctionParser::removeFromTypedValues", "will NOT remove " << (*s)->name());
	  
      }
      // remove all (or none if(symbolsToRemove.empty()) )
      for(typed_value_list_t::const_iterator s = symbolsToRemove.begin(); s != symbolsToRemove.end(); ++s)
	usedSymbols.remove(*s);
    }
  }

}
