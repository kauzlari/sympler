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




/*
Parts of the code in this file are based on code from E. Rudnyi's tdlib http://evgenii.rudnyi.ru/soft/tdlib00+/
*/


#include <math.h>
#include <string>
#include <cstdio>

#include "calc.h"
#include "general.h"

using namespace std;

calculator::map_entity *calculator::global;
calculator::map_key *calculator::delim_name;
bool calculator::qty = false;
calculator::Destruct calculator::clean;



double step_function(double x) {
  if (x > 0.0)
    return 1.0;
  else
    return 0.0;
}


double pow2(double x) {
  return x*x;
}


void calculator::init()
{
    if (!qty)
    {
        qty = true;
        global = new map_entity;
        typedef map_entity::value_type value;
        (*global).insert(value("acos", calculator_entity(acos)));
        (*global).insert(value("abs", calculator_entity(fabs)));
        (*global).insert(value("asin", calculator_entity(asin)));
        (*global).insert(value("atan", calculator_entity(atan)));
        (*global).insert(value("cos", calculator_entity(cos)));
        (*global).insert(value("exp", calculator_entity(exp)));
        (*global).insert(value("log10", calculator_entity(log10)));
        (*global).insert(value("log", calculator_entity(log)));
        (*global).insert(value("sin", calculator_entity(sin)));
        (*global).insert(value("sqrt", calculator_entity(sqrt)));
        (*global).insert(value("tan", calculator_entity(tan)));
        (*global).insert(value("R", calculator_entity(global::R)));
        (*global).insert(value("pi", calculator_entity(M_PI)));
        (*global).insert(value("step", calculator_entity(step_function)));
        (*global).insert(value("pow2", calculator_entity(pow2)));

        delim_name = new map_key;
        (*delim_name)['+'] = PLUS;
        (*delim_name)['-'] = MINUS;
        (*delim_name)['*'] = MUL;
        (*delim_name)['/'] = DIV;
        (*delim_name)['^'] = POW;
        (*delim_name)['('] = LP;
        (*delim_name)[')'] = RP;
        (*delim_name)[','] = COMMA;
    }
}

void calculator::get_token()
{
    *in_ptr >> ws;
    char t;
map_key::iterator i;
    in_ptr->get(t);
    if (!*in_ptr)
    {
        curr_tok.key = END;
        return;
    }
    if ((i = delim_name->find(t)) != delim_name->end())
    {
        text += t;
        curr_tok.key = (*i).second;
        return;
    }
    else if (isdigit(t) || t == '.')
    {
        in_ptr->putback(t);
        *in_ptr >> curr_tok.x;
        curr_tok.key = DOUBLE;
        char buf[30];
        sprintf(buf, "%f", curr_tok.x);
        text += buf;
        //    text += gcvt(curr_tok.x, 7, buf);
        return;
    }
    else if (t == ';')
    {
        curr_tok.key = END;
        return;
    }
    else if (isalnum(t) || t == '_' || t == '[' || t == ']' || t == '{' || t == '}')
    {
        string id(1, t);
        while (in_ptr->get(t), *in_ptr && (isalnum(t) || t == '_' || t == '[' || t == ']'
                || t == '{' || t == '}'
                || t == '.'))
          id += t;
        in_ptr->putback(t);
        map_entity::iterator i;
        text += id;
        if ((i = global->find(id)) != global->end())
            curr_tok = (*i).second;
        else {
            curr_tok.tag = id;
            AnalyzeId(id);
        }
    }
    else
      {
	string dummy(&t);
	throw gError("calculator::get_token", "unknown token \"" + dummy + "\"");
      }
}

void calculator::expr(bool must)
{
    calculator_entity tmp;
    term(must);
    for(;;)
        switch ((tmp = curr_tok).key) {
            case PLUS: case MINUS:
                get_token();
                term(true);
                tokens.push_back(tmp);
                break;
            case END: case RP:
                return;
            default:
                throw gError("calculator::expr", "NOTERM");
        }
}

void calculator::term(bool must)
{
    calculator_entity tmp;
    fact(must);
    for(;;)
        switch ((tmp = curr_tok).key) {
            case MUL: case DIV:
                get_token();
                fact(true);
                tokens.push_back(tmp);
                break;
            case END: case PLUS: case MINUS: case RP:
                return;
            default:
                throw gError("calculator::term", "NOFACT");
        }
}

void calculator::fact(bool must)
{
    calculator_entity tmp;
    prim(must);
    for(;;)
        switch ((tmp = curr_tok).key) {
            case POW:
                get_token();
                prim(true);
                tokens.push_back(tmp);
                break;
            case END: case PLUS: case MINUS: case MUL: case DIV: case RP:
                return;
            default:
                throw gError("calculator::fact", "NOPRIM");
        }
}

void calculator::prim(bool must)
{
    calculator_entity tmp;
    switch (curr_tok.key)
    {
        case DOUBLE:
        case DOUBLE_PTR:
        case DOUBLE_PTR_PTR:
            tokens.push_back(curr_tok);
            get_token();
            return;
        case LP:
            get_token();
            expr(true);
            if (curr_tok.key != RP) throw gError("calculator::prim",  "NORP");
                get_token();
            return;
        case FUNC1:
            tmp = curr_tok;
            get_token();
            if (curr_tok.key != LP) throw gError("calculator::prim", "NOLP");
                get_token();
            expr(true);
            if (curr_tok.key != RP) throw gError("calculator::prim", "NORP");
                get_token();
            tokens.push_back(tmp);
            return;
        case FUNC2:
            tmp = curr_tok;
            get_token();
            if (curr_tok.key != LP) throw gError("calculator::prim", "NOLP");
                get_token();
            expr(true);
            if (curr_tok.key != COMMA) throw gError("calculator::prim", "NOCOMMA");
                get_token();
            expr(true);
            if (curr_tok.key != RP) throw gError("calculator::prim", "NORP");
                get_token();
            tokens.push_back(tmp);
            return;
        case FUNC3:
            tmp = curr_tok;
            get_token();
            if (curr_tok.key != LP) throw gError("calculator::prim", "NOLP");
                get_token();
            expr(true);
            if (curr_tok.key != COMMA) throw gError("calculator::prim", "NOCOMMA");
                get_token();
            expr(true);
            if (curr_tok.key != COMMA) throw gError("calculator::prim", "NOCOMMA");
                get_token();
            expr(true);
            if (curr_tok.key != RP) throw gError("calculator::prim", "NORP");
                get_token();
            tokens.push_back(tmp);
            return;
        case PLUS:
            get_token();
            prim(true);
            return;
        case MINUS:
            tmp.key = UMIN;
            get_token();
            prim(true);
            tokens.push_back(tmp);
            return;
        case RP:
            throw gError("calculator::prim", "NOLP");
        case END:
            if (must) throw gError("calculator::prim", "NOPRIM");
            return;
        default:
            throw gError("calculator::prim", "NOPRIM");
    }
}

double calculator::eval() const
{
    if (tokens.empty()) return 0.;
    double arg2, arg3;
    vec_entity::const_iterator i;
    for (i = tokens.begin(); i != tokens.end(); i++)
    {
        switch ((*i).key)
        {
            case DOUBLE:
                stack.push_back((*i).x);
                break;
            case DOUBLE_PTR:
                stack.push_back(*(*i).x_ptr);
                break;
            case DOUBLE_PTR_PTR:
                stack.push_back(**(*i).x_ptr_ptr);
                break;
            case FUNC1:
                stack.back() = (*(*i).f1)(stack.back());
                break;
            case FUNC2:
                arg2 = stack.back();
                stack.pop_back();
                stack.back() = (*(*i).f2)(stack.back(), arg2);
                break;
            case FUNC3:
                arg3 = stack.back();
                stack.pop_back();
                arg2 = stack.back();
                stack.pop_back();
                stack.back() = (*(*i).f3)(stack.back(), arg2, arg3);
                break;
            case PLUS:
                arg2 = stack.back();
                stack.pop_back();
                stack.back() += arg2;
                break;
            case MINUS:
                arg2 = stack.back();
                stack.pop_back();
                stack.back() -= arg2;
                break;
            case MUL:
                arg2 = stack.back();
                stack.pop_back();
                stack.back() *= arg2;
                break;
            case DIV:
                arg2 = stack.back();
                stack.pop_back();
                stack.back() /= arg2;
                break;
            case POW:
                arg2 = stack.back();
                stack.pop_back();
                stack.back() = pow(stack.back(), arg2);
                break;
            case UMIN:
                stack.back() = -stack.back();
                break;
            default:
                throw gError("calculator::est", "run-time WRTOK");
        }
    }
    if (stack.size() != 1)
        throw gError("calculator::est", "run-time WRONG_TOKENS");
    arg2 = stack.back();
    stack.pop_back();
    return arg2;
}


bool calculator::isConstant() const
{
  bool is_constant = true;

  if (!tokens.empty()) {
    vec_entity::const_iterator i;
    for (i = tokens.begin(); i != tokens.end() && is_constant; i++)
    {
      if (i->key == DOUBLE_PTR || i->key == DOUBLE_PTR_PTR)
        is_constant = false;
    }
  }

  return is_constant;
}


set<string> calculator::enumVariables() const
{
  set<string> vars;

  if (!tokens.empty()) {
    vec_entity::const_iterator i;
    for (i = tokens.begin(); i != tokens.end(); i++)
    {
      if (i->key == DOUBLE_PTR || i->key == DOUBLE_PTR_PTR)
        vars.insert(i->tag);
    }
  }

  return vars;
}



//---- class fn ----

void fn::AnalyzeId(const string &id)
{
    map<string, double>::iterator i;

    if (id == "x") {
        curr_tok.key = DOUBLE_PTR;
        curr_tok.x_ptr = &m_x;
    } else if ((i = m_constants.find(id)) != m_constants.end()) {
        curr_tok.key = DOUBLE;
        curr_tok.x = i->second;
    } else
        throw gError("fn::AnalyzeId", "Id undefined. ");
}

