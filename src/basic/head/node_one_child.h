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



#ifndef __NODE_ONE_CHILD_H
#define __NODE_ONE_CHILD_H

#include "node.h"

#include "general.h"

/*!
 * The base class of all object that can have exactly one child.
 */
class NodeOneChild: public Node
{
protected:
  /*!
   * The pointer to the child
   */
  Node *m_child;

  /*!
   * Request instantiation of a child named \a name. This is called
   * during the XML read procedure, i.e., the \a read method
   * @param name The name of the node
   */
  virtual Node *instantiateChild(const string &name) = 0;

public:
  /*!
   * Constructor
   */
  NodeOneChild();

  /*!
   * Constructor
   * @param parent Pointer to the parent node
   */
  NodeOneChild(Node *parent);

  /*!
   * Destructor
   */
  virtual ~NodeOneChild();

  /*!
   * \a read reads all properties from the XML file
   * @param xmln The corresponding xmln (uses libXML2)
   */
  virtual void read(const xmlNode *xmln);

  /*!
   * \a write writes all properties to an XML file
   * @param s Output file stream
   * @param shift Number of spaces to be introduced in the beginning of the line
   */
  virtual ostream &write(ostream &s, int shift = 0);

  /*!
   * \a setup initializes all the variables. This ensures object
   * creation of all entities in the XML file is complete when setup is being called.
   * It also invokes setup for all child objects.
   */
  virtual void setup();
};

#endif
