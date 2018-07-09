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



#include "node_one_child.h"

//---- Construtors/Destructor ----

NodeOneChild::NodeOneChild(): Node(), m_child(NULL)
{
}


NodeOneChild::NodeOneChild(Node *parent): Node(parent), m_child(NULL)
{
}


NodeOneChild::~NodeOneChild()
{
  delete m_child;
}



//---- Methods ----

/*
void NodeOneChild::read(const SGML &e)
{
  Node::read(e);

  parser p1(e);
  SGML e2;

  while (!p1.eof()) {
    p1.GetSGML(e2);

    if (m_child)
      throw gError("NodeOneChild::read", "Only one child is allowed for this node.");
    else {
      m_child = instantiateChild(e2.name);
      if (m_child)
        m_child->read(e2);
    }
  }    
}
*/


void NodeOneChild::read(const xmlNode *xmln)
{
  Node::read(xmln);

  for (const xmlNode *cur_node = xmln->children; cur_node; cur_node = cur_node->next) {
    if (cur_node->type == XML_ELEMENT_NODE) {
      if (m_child)
        throw gError("NodeOneChild::read", "For " + className() + ": Only one child is allowed for this node.");

      m_child = instantiateChild((const char*) cur_node->name);
      if (m_child)
        m_child->read(cur_node);
    }
  }
}


ostream &NodeOneChild::write(ostream &s, int shift)
{
  m_properties.toXML_begin(s, shift);
    
  if (m_child)
    m_child->write(s, shift+1);

  writeMoreBody(s, shift);
    
  m_properties.toXML_end(s, shift);
    
  return s;
}


void NodeOneChild::setup()
{
  Node::setup();

  if (m_child) {
    MSG_DEBUG("NodeOneChild::setup", "Running setup for: " << m_child->className());
    m_child->setup();
  }
}
