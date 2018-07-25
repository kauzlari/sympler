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



#include <algorithm>
#include <functional>

#include "node_many_children.h"

#include "general.h"

// #include "valgrind/memcheck.h"

//---- Construtors/Destructor ----

NodeManyChildren::NodeManyChildren(): Node()
{
}


NodeManyChildren::NodeManyChildren(Node *parent): Node(parent)
{
}


NodeManyChildren::~NodeManyChildren()
{
  for_each(m_children.begin(), m_children.end(), delete_node);
}



//---- Methods ----


struct find_child: public unary_function<const string&, bool> {
  string m_name;
  find_child(const string &name): m_name(name) { }
  bool operator()(Node *c) {
//     MSG_DEBUG("find_child", "c->className() = '" << c->className() << "', name = '" << m_name << "', c->name() = '" << c->name() << "'");

    return c->/*name*/className() == m_name;
  }
};


Node *NodeManyChildren::findChild(const string &name)
{
  list<Node*>::iterator n = find_if(m_children.begin(), m_children.end(), find_child(name));

  if (n == m_children.end())
    return NULL;
  else
    return *n;
}

/*
void NodeManyChildren::read(const SGML &e)
{
  Node::read(e);

  parser p1(e);
  SGML e2;

  while (!p1.eof()) {
    Node *child;
        
    p1.GetSGML(e2);

    child = instantiateChild(e2.name, e2);
    if (child) {
      child->read(e2);

      m_children.push_back(child);
    }
  }
}
*/


void NodeManyChildren::read(const xmlNode *xmln)
{
  Node::read(xmln);

  for (const xmlNode *cur_node = xmln->children; cur_node; cur_node = cur_node->next) {
    if (cur_node->type == XML_ELEMENT_NODE) {
      Node *child;

      child = instantiateChildWithXML((const char*) cur_node->name, cur_node);
      if (child) {
        child->read(cur_node);

        m_children.push_back(child);
      }
    }
  }
}


ostream &NodeManyChildren::write(ostream &s, int shift)
{
  m_properties.toXML_begin(s, shift);
    
  for (list<Node*>::iterator i = m_children.begin(); i != m_children.end(); i++)
    (*i)->write(s, shift+1);

  writeMoreBody(s, shift);
    
  m_properties.toXML_end(s, shift);
    
  return s;
}


void NodeManyChildren::setup()
{
  Node::setup();

//    MSG_DEBUG("NodeManyChildren::setup","for " << this->className() << ": VALGRINDSTART:");
//   MSG_DEBUG("NodeManyChildren::setup", "Valgrindcheck m_children: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_children));

  FOR_EACH
    (list<Node*>,
     m_children,
     MSG_DEBUG("NodeManyChildren::setup", "Running setup for: " << (*__iFE)->className());
     (*__iFE)->setup();
    );
}

