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



#ifndef __VARIANT_H
#define __VARIANT_H

#include <vector>

using namespace std;

#include <assert.h>

#include "general.h"
#include "geometric_primitives.h"


typedef vector<double> vector_double_t;
typedef vector<string> vector_string_t;


/*!
 * The Variant class contains both, values and type
 * information which are determined during runtime.
 */
class Variant
{
 public:
  enum variant_type_t {
    SCALAR = 0,
    VECTOR = 1,
    TENSOR = 2,
    SCALAR_STRING = 3,
    VECTOR_STRING = 4,
    TENSOR_STRING = 5,
    EO_TYPED_VALUE = 6
  };

 protected:
  variant_type_t m_type_id;

  vector_double_t m_doubles;
  vector_string_t m_strings;

 public:
  Variant(variant_type_t t);
  Variant(const Variant &copy);
  ~Variant();

  Variant &operator=(const Variant &copy);

  variant_type_t typeId() const {
    return m_type_id;
  }


  /* double types */
  double &scalar() {
    assert(m_type_id == SCALAR);

    return m_doubles[0];
  }
  double scalar() const {
    assert(m_type_id == SCALAR);

    return m_doubles[0];
  }

  double &vector(int i) {
    assert(m_type_id == VECTOR);

    return m_doubles[i];
  }
  double vector(int i) const {
    assert(m_type_id == VECTOR);

    return m_doubles[i];
  }

  double &tensor(int i, int j) {
    assert(m_type_id == TENSOR);

    return m_doubles[i*SPACE_DIMS+j];
  }
  double tensor(int i, int j) const {
    assert(m_type_id == TENSOR);

    return m_doubles[i*SPACE_DIMS+j];
  }

  vector_double_t &doubles() {
    return m_doubles;
  }

  /* string types */
  string &scalarString() {
    assert(m_type_id == SCALAR_STRING);

    return m_strings[0];
  }
  string scalarString() const {
    assert(m_type_id == SCALAR_STRING);

    return m_strings[0];
  }

  string &vectorString(int i) {
    assert(m_type_id == VECTOR_STRING);

    return m_strings[i];
  }
  string vectorString(int i) const {
    assert(m_type_id == VECTOR_STRING);

    return m_strings[i];
  }

  string &tensorString(int i, int j) {
    assert(m_type_id == TENSOR_STRING);

    return m_strings[i*SPACE_DIMS+j];
  }
  string tensorString(int i, int j) const {
    assert(m_type_id == TENSOR_STRING);

    return m_strings[i*SPACE_DIMS+j];
  }

  vector_string_t &strings() {
    return m_strings;
  }

  /* Simple checks which bail out if they fail */
  variant_type_t sameType(const Variant &b) const; /*{
    if (m_type_id != b.m_type_id) {
      throw gError
	("Variant::sameType",
	 "Type mismatch.");
    }

    return m_type_id;
  }*/
};



#endif
