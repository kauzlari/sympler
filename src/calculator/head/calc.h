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


#ifndef __CALC_H
#define __CALC_H

//#include "SGML.h"

#include "general.h"

using namespace std;

inline string RemoveSpaces(const string &str)
{
	string out;
	for (string::const_iterator i = str.begin(); i != str.end(); ++i)
		if (!isspace(*i))
			out.append(1, *i);
	return out;
}

class calculator
{
protected:
  enum calculator_key {DOUBLE, DOUBLE_PTR, DOUBLE_PTR_PTR, FUNC1, FUNC2, FUNC3,
                       PLUS, MINUS, MUL, DIV, POW, LP, RP, COMMA, UMIN,
                       END};
  class calculator_entity;
  friend class calculator::calculator_entity;
  class calculator_entity
  {
  public:
    calculator_key key;
    union
    {
      double x;
      double *x_ptr;
      double **x_ptr_ptr;
      double (*f1)(double);
      double (*f2)(double, double);
      double (*f3)(double, double, double);
    };
    string tag;

    calculator_entity() {key = END;}
    calculator_entity(double x_) {key = DOUBLE; x = x_;}
    calculator_entity(double *x_ptr_) {key = DOUBLE_PTR; x_ptr = x_ptr_;}
    calculator_entity(double **x_ptr_ptr_)
      {key = DOUBLE_PTR_PTR; x_ptr_ptr = x_ptr_ptr_;}
    calculator_entity(double (*f1_)(double)) {key = FUNC1; f1 = f1_;}
    calculator_entity(double (*f2_)(double, double))
      {key = FUNC2; f2 = f2_;}
    calculator_entity(double (*f3_)(double, double, double))
      {key = FUNC3; f3 = f3_;}
  };

  typedef vector<calculator_entity> vec_entity;
  typedef map<string, calculator_entity, less<string> > map_entity;
  typedef map<char, calculator_key, less<char> > map_key;

  static map_entity *global;
  static map_key *delim_name;
  static bool qty;

  vec_entity tokens;
  string text;
  mutable vec_double stack;

  calculator_entity curr_tok;
  istream* in_ptr;

  virtual void init();

  virtual void get_token();
  virtual void expr(bool must);
  virtual void term(bool must);
  virtual void fact(bool must);
  virtual void prim(bool must);

  virtual void AnalyzeId(const string &id) = 0;
//  virtual string AnalyzeSGML(const SGML &e);

private:
  class Destruct;
  friend class Destruct;
  class Destruct
  {
  public:
    ~Destruct()
    {
        if (calculator::qty)
        {
            calculator::qty = false;
            delete calculator::global;
            delete calculator::delim_name;
        }
    }
  };
  static Destruct clean;

public:

  calculator() {init();}
  calculator(const calculator& old) {init(); FromString(old.text);}
  calculator(const string &old) {init(); FromString(old);}
  calculator& operator=(const calculator& old)
  {
		if (this != &old) 
			FromString(old.text); 
		return *this;
	}
  virtual ~calculator() {}

  virtual void FromString(const string& old)
  {
    istringstream in(old);
    read(in);
  }
  virtual void clear()
  {
    tokens.clear();
    text.erase();
  }

  virtual istream& read(istream& in)
  {
    clear();
    in_ptr = &in;
    get_token();
    expr(false);
    stack.reserve(10);
    return in;
  }
  virtual ostream& write(ostream& out) const
    {return out << text;}
  virtual const string& expression() const
    {return text;}
  virtual bool empty()
    {return text.empty();}

  virtual double eval() const;

  virtual bool isConstant() const;
  virtual set<string> enumVariables() const;
};

OVERLOAD_STREAMS(calculator)


class fn: public calculator
{
protected:
    double m_x;
    map<string, double> m_constants;

    virtual void AnalyzeId(const string &id);

public:
    fn(): calculator(), m_x(0.0) { }
    fn(const fn& old): calculator() {
        FromString(old.text);
    }
    fn(const string& old): calculator() {
        FromString(old);
    }

    virtual double &constant(const string &name) {
        return m_constants[name];
    }

    virtual double &x() {
        return m_x;
    }

    virtual double value(double x) {
        m_x = x;
        return eval();
    }
};


#endif

