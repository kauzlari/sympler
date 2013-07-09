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



#ifndef __FUNCTION_NODE_H
#define __FUNCTION_NODE_H

#include <string>

#include "variant.h"
#include "smart_enum.h"

using namespace std;

#define C_MAX_PRIORITY 10

/*!
 * A FunctionNode is the base class for the nodes in the 
 * parse tree.
 */
class FunctionNode
{
 public:
  enum node_type_t {
    OTHER = 0,
    UNARY_FUNCTION = 1,
    BINARY_FUNCTION = 2,
    UNARY_OPERATOR = 3,
    BINARY_OPERATOR = 4,
    EO_NODE_TYPE = 5
  };

  /*!
   * The name/symbol of this operation/variable/etc.
   */
  string m_name;

 public:
  /*!
   * The constructor
   */
  FunctionNode(string name);

  /*!
   * The destructor 
   */
  virtual ~FunctionNode();

  /*!
   * Return the name of this node
   */
  virtual string name() const {
//     MSG_DEBUG("FunctionNode::name", "name = " << m_name);
    return m_name;
  }

  /*!
   * Return the value
   */
  virtual Variant value() const = 0;

  /*!
   * Output C code for external compilation
   */
  virtual Variant toC() const = 0;
};


/*!
 * The factory producing function nodes.
 */
class FunctionNode_Factory: public SmartEnum<FunctionNode_Factory>
{
protected:
  FunctionNode::node_type_t m_node_type;

  static list<FunctionNode_Factory*> *s_nodes_by_priority;

  static list<FunctionNode_Factory*> &priorityF(int priority) {
    if (!s_nodes_by_priority)
      s_nodes_by_priority = new list<FunctionNode_Factory*>[C_MAX_PRIORITY];

    return s_nodes_by_priority[priority];
  }

  FunctionNode_Factory(FunctionNode::node_type_t t, int priority, const string &name)
    : SmartEnum<FunctionNode_Factory>(name), m_node_type(t) {
    priorityF(priority).push_back((FunctionNode_Factory*) this);
  }

public:
  virtual FunctionNode *instantiate() const = 0;

  virtual FunctionNode::node_type_t nodeType() {
    return m_node_type;
  }

  static const list<FunctionNode_Factory*>& byPriority(int priority) {
    return s_nodes_by_priority[priority];
  }
};


template <class T>
class FunctionNode_Register: public FunctionNode_Factory
{
 public:
  FunctionNode_Register(FunctionNode::node_type_t t, int priority, const string &symbol)
    : FunctionNode_Factory(t, priority, symbol) {
  }
    
  virtual FunctionNode *instantiate() const {
    return new T();
  }
};



/* Macros */

#define FOR_EACH_DOUBLE(where, varianta, code)				\
  {									\
    Variant::variant_type_t t;						\
    Variant vca = varianta;						\
									\
    t = vca.typeId();							\
									\
    if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) { \
      Variant r(t);							\
      vector_double_t &vr = r.doubles();				\
      vector_double_t va = vca.doubles();				\
									\
      for (size_t i = 0; i < va.size(); ++i) {				\
	code;								\
      }									\
									\
      return r;								\
    } else								\
      throw gError							\
	(#where,							\
	 "Can't handle Variant type id " + ObjToString(t) + ".");	\
  }


#define FOR_EACH_STRING(where, varianta, code)				\
  {									\
    Variant::variant_type_t t;						\
    Variant vca = varianta;						\
									\
    t = vca.typeId();							\
									\
    if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	\
	t == Variant::TENSOR_STRING) {					\
      Variant r(t);							\
      vector_string_t &vr = r.strings();				\
      vector_string_t va = vca.strings();				\
									\
      for (size_t i = 0; i < va.size(); ++i) {				\
	code;								\
      }									\
									\
      return r;								\
    } else								\
      throw gError							\
	(#where,							\
	 "Can't handle Variant type id " + ObjToString(t) + ".");	\
  }

#endif
