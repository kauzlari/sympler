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




/*
Parts of the code in this file are based on code from E. Rudnyi's tdlib http://evgenii.rudnyi.ru/soft/tdlib00+/
*/


#ifndef __GENERAL_H
#define __GENERAL_H

#include <map>
#include <set>
#include <list>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

#include <assert.h>
#include <pthread.h>

#include "misc.h"
#include "consts.h"

#include<typeinfo>
#include<algorithm>
#include<string.h>

#include <cstdio>

#ifdef HAVE_TNT_TNT_H
#define WITH_ARRAY_TYPES
#endif

using namespace std;


//---- Macros ----

#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif

#define MSG_DEBUG(where, what)  cout << "{" << where << "} " << what << endl
#define MSG_INFO(where, what)  MSG_DEBUG(where, what)
#define __LINE_STRING__ ObjToString(__LINE__)
#define SHORTFILE (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define FILE_INFO (std::string("(") + SHORTFILE + ":" + __LINE_STRING__ + ")")
#define FILEFUNC_INFO (std::string("") + __func__ + " (" + SHORTFILE + ":" + __LINE_STRING__ + ")")
#ifdef __PRETTY_FUNC__
#define PRETTYFUNC_INFO (std::string("") + __PRETTY_FUNC__ + " (" + SHORTFILE + ":" + __LINE_STRING__ + ")")
#else
#define PRETTYFUNC_INFO FILEFUNC_INFO
#endif

//---- Typedefs ----

typedef set<int> group_t;
typedef vector<point_t> forces_t;

// this is a definition for iterators over particles; if the array type of my
// particle array changes in phase.h, this here has also to be changed
//typedef list<int>::iterator PIter;
typedef vector<double> vec_double;
typedef vec_double::iterator vec_double_i;
typedef vec_double::const_iterator vec_double_ci;
typedef vector<double*> vec_double_ptr;
typedef vec_double_ptr::iterator vec_double_ptr_i;
typedef vec_double_ptr::const_iterator vec_double_ptr_ci;
typedef vector<double**> vec_double_2ptr;
typedef vec_double_2ptr::iterator vec_double_2ptr_i;
typedef vec_double_2ptr::const_iterator vec_double_2ptr_ci;
typedef vector<string> vec_string;
typedef vec_string::iterator vec_string_i;
typedef vec_string::const_iterator vec_string_ci;
typedef vector<int> vec_int;
typedef vec_int::iterator vec_int_i;
typedef vec_int::const_iterator vec_int_ci;
typedef map<string, string, less<string> > map_string_string;
typedef map_string_string::iterator map_string_string_i;
typedef map_string_string::const_iterator map_string_string_ci;


#define index index_


//the idea of this class is taken from gnussl
class gError
{
protected:
  string m_where;		// where did the error occur?
  string m_msg;		// the message
  int m_exitValue;  // exit value

public:
	enum exitValues {
		DEFAULT = 1,
		PARTICLEFLEWTOOFAR = 2
	};

  gError() {}
  gError(const char *msg) : m_where("Unknown"), m_msg(string("ERROR: ") + msg), m_exitValue(DEFAULT) {}
  gError(const string msg) : m_where("Unknown"), m_msg(string("ERROR: ") + msg), m_exitValue(DEFAULT) {}
  gError(const char *where, const char *msg) : m_where(where), m_msg(string("ERROR: ") + msg), m_exitValue(DEFAULT) { }
  gError(const string where, const string msg) : m_where(where), m_msg(string("ERROR: ") + msg), m_exitValue(DEFAULT) { }
  gError(const char *msg, const int exitValue) : m_where("Unknown"), m_msg(string("ERROR: ") + msg), m_exitValue(exitValue) {}
  gError(const string msg, const int exitValue) : m_where("Unknown"), m_msg(string("ERROR: ") + msg), m_exitValue(exitValue) {}
  gError(const char *where, const char *msg, int exitValue) : m_where(where), m_msg(string("ERROR: ") + msg), m_exitValue(exitValue) { }
  gError(const string where, const string msg, int exitValue) : m_where(where), m_msg(string("ERROR: ") + msg), m_exitValue(exitValue) { }

  virtual ~gError() { }

  virtual string where() {
    return m_where;
  }

  virtual string msg() {
    return m_msg;
  }

  virtual string message() {         // the error message
    if (m_where == "Unknown")
      return m_msg;

    return "[" + m_where + "] " + m_msg;
  }

  virtual int exitValue() {
	  return m_exitValue;
  }
};


#define OVERLOAD_STREAMS(T)                     \
inline istream& operator>>(istream &in, T &to)  \
{                                               \
  in >> ws;                                     \
  return to.read(in);                           \
}                                               \
                                                            \
inline ostream& operator<<(ostream &out, const T &old)      \
{                                                           \
  return old.write(out);                                    \
}


class PutTab
{
public:
  size_t n;
  PutTab(size_t n) : n(n) {}
};

inline ostream& operator<<(ostream &out, const PutTab &old)
{
  for (size_t i = 0; i < old.n; ++i)
    out << "  ";
  return out;
}

template <class T>
string ObjToString(const T &x)
{
  ostringstream out;
  out << x;
  string res(out.str());
  return res;
}


class global
{
public:
  static double R;     // R = 8.31441 by default
  static size_t n_threads;
#ifdef _MPI
  static size_t thread_id;
  static size_t n_procs;
#endif
};


/*!
 * Helper function checking if string \a s is within list of strings
 * \a pipeStringList . The string \a pipeStringList is assumed to have the form
 * "string1|string2|...|stringN".
 *
 */
bool g_stringIsInPipeList(string s, string pipeStringList);


/*! make_filename takes the filename 's' and appends the number 'n'
   in front of the suffix. */
string make_filename(string s, int n);

#endif

