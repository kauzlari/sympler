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





/* The code in this file is an adapted version of a template class 
 * for quasi-enum classes originally implemented in the LibBase 
 * Source Code (for Copyright see below) which was obtained with kind 
 * permission from Jan G. Korvink.
 */

/* LibBase Source Code
   Copyright by Markus Emmenegger, 
   Stefano Taschini, Jan G. Korvink, 1996,1997,1998
   Physical Electronics Laboratory, ETH-Hoenggerberg, CH-8093 Zuerich
   Email:   solidis@iqe.phys.ethz.ch
   All rights reserved
   template class for quasi-enum classes, meant to replace the enums
*/



/* Some compiler/linker cannot handle SmartEnum correctly. In this case
   try to define ALTERNATIVE_SMART_ENUMS and see if it works. */


#ifndef __SMART_ENUM_H
#define __SMART_ENUM_H 

/* We choose our SmartEnum method for now because it works on AIX and links
   properly. */
#define ALTERNATIVE_SMART_ENUM

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include "general.h"

using namespace std;

// This class is intended to replace enums
// For instance instead of
// enum myType {aaa,bbb,ccc};
// one can use
/*
 class myType: public LBSmartEnum<class myType>
 {
public:
     static const myType aaa;
     static const myType bbb;
     static const myType ccc;
 };
 */


#ifdef ALTERNATIVE_SMART_ENUM
//#define mapF()  s_map
//#define vectorF()  s_list
#endif


/*!
 * Class database with appropriate object factories
 */
template <class T>
class SmartEnum
{
protected:
  /*!
   * Name of this class
   */
  string fName;
  
  /*!
   * Index of this class
   */
  int fOrdinal;

  /*!
   * Description of this database
   */
  static string s_description;

#ifdef ALTERNATIVE_SMART_ENUM
  static map<string, const T*, less<string> > *s_map;
  static vector<const T*> *s_list;

  static map<string, const T*, less<string> >& mapF() {
    if (!s_map) {
//      MSG_DEBUG("SmartEnum::mapF", "Creating map");

      s_map = new map<string, const T*, less<string> >;
    }

    return *s_map;
  }
    
  static vector<const T*>& vectorF() {
    if (!s_list) {
//      MSG_DEBUG("SmartEnum::vectorF", "Creating vector");

      s_list = new vector<const T*>;
    }

    return *s_list;
  }
#else
  static map<string, const T*, less<string> >& mapF() {
    static map<string, const T*, less<string> > aMap;
    return aMap;
  }
    
  static vector<const T*>& vectorF() {
    static vector<const T*> aList;
    return aList;
  }
#endif

  /*!
   * Add factory for class \a name to this database
   * @param name Name of the new class
   */
  SmartEnum(const string& name) {
    fName = name;
    fOrdinal = cardinality();

    if (!mapF().insert(make_pair<const string, const T*>(name, (T*) this)).second) {
      throw gError("SmartEnum::SmartEnum", "Unable to register object " + name + " .");
    }

    vectorF().push_back((T*) this);
  }

  /*!
   * Destructor
   */
  virtual ~SmartEnum() {
  }

public:
  /*!
   * Return the name of this class
   */
  const string& name() const {
    return fName;
  }

  /*!
   * Return the index of this class
   */
  int ordinal() const {
    return fOrdinal;
  }
    
  /*!
   * Return the factory for class \a name
   * @param name Class name
   */
  static const T& byName(const string& name) {
    typename map<string, const T*, less<string> >::iterator p = mapF().find(name);

    if(p == mapF().end()) {
      throw gError("SmartEnum::byName", "Object " + name + " not found in database.");
    }

    return *(p->second);
  }

  /*!
   * Does a class of name \a name exist in this database?
   * @param name Class name
   */
  static bool exists(const string &name) {
    typename map<string, const T*, less<string> >::iterator p = mapF().find(name);

    return p != mapF().end();
  }

  /*!
   * Return the number of entries in this class database
   */
  static int cardinality() {
    return vectorF().size();
  }
    
  /*!
   * Return factory by index
   * @param ordinal Factory index
   */
  static const T& byOrdinal(int ordinal) {
    return *vectorF()[ordinal];
  }

  /*!
   * Compare the two factories
   * @param op The other factory
   */
  bool operator==(const SmartEnum<T>& op) const {
    return this == &op;
  }
    
  /*!
   * Compare the two factories
   * @param op The other factory
   */
  bool operator!=(const SmartEnum<T>& op) const {
    return this != &op;
  }

  /*!
   * Return the description of this database
   */
  static string description() {
    return s_description;
  }    
};


#ifdef ALTERNATIVE_SMART_ENUM
//#undef mapF
//#undef vectorF

#define REGISTER_SMART_ENUM(factory, description)                       \
  template<> string SmartEnum<factory>::s_description(description);                \
  template<> map<string, const factory*, less<string> > *SmartEnum<factory>::s_map = NULL; \
  template<> vector<const factory*> *SmartEnum<factory>::s_list = NULL


#else

#define REGISTER_SMART_ENUM(factory, description)       \
template<> string SmartEnum<factory>::s_description(description)

#endif

#endif
