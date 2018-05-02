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


#include "variant.h"

#include "general.h"


Variant::Variant(variant_type_t t): m_type_id(t)
{
  switch (m_type_id) {
  case SCALAR:
    m_doubles.resize(1, 0);
    break;
  case VECTOR:
    m_doubles.resize(SPACE_DIMS, 0);
    break;
  case TENSOR:
    m_doubles.resize(SPACE_DIMS*SPACE_DIMS, 0);
    break;
  case SCALAR_STRING:
    m_strings.resize(1, "");
    break;
  case VECTOR_STRING:
    m_strings.resize(SPACE_DIMS, "");
    break;
  case TENSOR_STRING:
    m_strings.resize(SPACE_DIMS*SPACE_DIMS, "");
    break;
  default:
    throw gError
      ("Variant::Variant",
       "(Internal error) Unknown variant type " + ObjToString(t) + ".");
  }
}


Variant::Variant(const Variant &copy): m_type_id(copy.m_type_id)
{
  m_doubles = copy.m_doubles;
  m_strings = copy.m_strings;
}


Variant::~Variant()
{
}


Variant &Variant::operator=(const Variant &copy)
{
  m_type_id = copy.m_type_id;
  m_doubles = copy.m_doubles;
  m_strings = copy.m_strings;

  return *this;
}

Variant::variant_type_t Variant::sameType(const Variant &b) const 
{
  if (m_type_id != b.m_type_id) {
    throw gError
        ("Variant::sameType",
          "Type mismatch.");
  }

  return m_type_id;
}
