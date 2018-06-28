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



#ifndef __NODE_H
#define __NODE_H

#include <iostream>
#include <libxml/tree.h>

#include "property_list.h"

/* General class structure:
   - Node can have a parent
   - NodeOneChild can have excatly one child (and of course a parent)
   - NodeManyChildren can have as many children as wished
   - Every Node has a property list
*/

// This is the changed Version

//---- Macros ----

#define INTPC(name, var, min, comm)  m_properties.addProperty(#name, PropertyList::INT, &var, new PLCIntGreater(min), comm); while(0)
#define DOUBLEPC(name, var, min, comm)  m_properties.addProperty(#name, PropertyList::DOUBLE, &var, new PLCDoubleGreater(min), comm); while(0)
#define STRINGPC(name, var, comm)  m_properties.addProperty(#name, PropertyList::STRING, &var, NULL, comm); while(0)
#define STRINGPCINF(name, var, comm)  m_properties.addProperty(#name+string("InputFile"), PropertyList::STRING, &var, NULL, comm); while(0)
#define STRINGPCOUF(name, var, comm)  m_properties.addProperty(#name+string("OutputFile"), PropertyList::STRING, &var, NULL, comm); while(0)
#define BOOLPC(name, var, comm)  m_properties.addProperty(#name, PropertyList::BOOLEAN, &var, NULL, comm); while(0)
#define FUNCTIONPAIRPC(name, var, comm)  m_properties.addProperty(#name, PropertyList::FUNCTIONPAIR, &var, NULL, comm); while(0)
#define FUNCTIONFIXEDPC(name, var, comm)  m_properties.addProperty(#name, PropertyList::FUNCTIONFIXED, &var, NULL, comm); while(0)
#define POINTPC(name, var, comm)  m_properties.addProperty(#name, PropertyList::POINT, &var, NULL, comm); while(0)


//---- Classes ----

/*!
 * The base class of all objects that can be autoloaded/configured
 * from an XML file into a tree structure.
 */
class Node
{
protected:
  /*!
   * The parent node (NULL if none)
   */
  Node *m_parent;
    
  /*!
   * The property list holding all properties that can be
   * autoloaded from an XML file.
   */
  PropertyList m_properties;

  /*!
   * Write the body part of the XML file.
   */
  virtual void writeMoreBody(ostream &s, int shift = 0);

public:
  /*!
   * Constructor
   */
  Node();

  /*!
   * Constructor
   * @param parent Pointer to the parent node
   */
  Node(Node *parent);

  /*!
   * Destructor
   */
  virtual ~Node();

  /*!
   * Return the pointer to the parent node
   */
  virtual Node *parent() {
    return m_parent;
  }

  /*!
   * Returns the \a PropertyList of this \a Node
   */
  const PropertyList& returnProperties() const {
    return m_properties;
  }
  
  /*!
   * String identifier for this object
   */
  const string &name() const {
    return m_properties.name();
  }

  /*!
   * String identifier for this class. Should
   * be identical to using RTTI
   */
  const string &className() const {
    return m_properties.className();
  }

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
   */
  virtual void setup();

  /*!
   * Additional setup routine
   */
  virtual void setupAfterParticleCreation()
  {
    
  }
  
  /*!
   * Possibility for the Node to precompute stuff before \a Symbol s and \a Force s start computing
   */
  virtual void precompute()
  {
    
  }
  
  /*!
   * Returns a help tree
   */
  virtual void help(HelpNode *node) {
    return m_properties.help(node);
  }
};


//---- Function -----

/* Basically intended to be used in for_each, because unfortunately
   C++ doesn't have lambda expressions. */
void delete_node(Node *n);

#endif
