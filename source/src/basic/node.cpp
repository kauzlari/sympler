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



#include "node.h"

// for MSG_DEBUG
#include "general.h"

// #include "valgrind/memcheck.h"

//---- Constructors/Destructors ----

Node::Node(): m_parent(NULL)
{
}


Node::Node(Node *parent): m_parent(parent)
{
}


Node::~Node()
{
}



//---- Methods ----

void Node::read(const xmlNode *xmln)
{
  m_properties.fromXML(xmln);
}


ostream &Node::write(ostream &s, int shift)
{
  m_properties.toXML_begin(s, shift);

  writeMoreBody(s, shift);

  m_properties.toXML_end(s, shift);

  return s;
}


void Node::writeMoreBody(ostream &s, int shift)
{
}


void Node::setup()
{
//   MSG_DEBUG("Node::setup","for " << this->className() << ": VALGRINDSTART:");
//   MSG_DEBUG("Node::setup", "Valgrindcheck m_parent: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_parent));
//   MSG_DEBUG("Node::setup", "Valgrindcheck m_properties: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_properties));
}

//---- Functions ----

void delete_node(Node *n)
{
    delete n;
}
