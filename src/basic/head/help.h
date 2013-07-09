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



#ifndef __HELP_H
#define __HELP_H

#include <list>
#include <string>
#include <ostream>

#include "node.h"
#include "smart_enum.h"

using namespace std;


/*!
 * A node of an entry in the help system
 */
class HelpNode
{
 protected:
  /*!
   * Parent node
   */
  HelpNode *m_parent;

  /*!
   * Child nodes
   */
  list<HelpNode*> m_children;

  /*!
   * Help text of this node
   */
  string m_text;

 public:
  /*!
   * Constructor
   */
  HelpNode();

  /*!
   * Construct help node with parent \a parent and help text \a text
   * @param parent Parent node
   * @param text Help text
   */
  HelpNode(HelpNode *parent, string text);

  /*!
   * Destructor
   */
  virtual ~HelpNode();

  /*!
   * Add a child to this help node
   */
  virtual void addChild(HelpNode *child);

  /*!
   * Return the text
   */
  const string &text() const {
    return m_text;
  }

  /*!
   * Return the children
   */
  const list<HelpNode*> &children() const {
    return m_children;
  }
};


/*!
 * Format the text \a str with \a indent spaces to the left
 * to it does not exceed the width \a width
 * @param str String to format
 * @param indent Number of spaces to put on the left
 * @param width Maximum width of the text block
 */
string block(string str, size_t indent, size_t width);

/*!
 * Produce help output from node \a root, write it to stream \a s with
 * initially \a start_depth spaces to the left.
 * @param s Output stream
 * @param root Root help node
 * @param start_depth Initial number of spaces to the left
 */
void HelpFormatScreen(ostream &s, const HelpNode &root, int start_depth);



/* --- */

template <class T>
class Help
{
 protected:
    virtual Node *instantiate(const T &factory) const = 0;

 public:
    Help() {
    }

    HelpNode *asHelpNode(HelpNode *root) const {
        for (int i = 0; i < T::cardinality(); i++) {
            Node *n;

            n = instantiate(T::byOrdinal(i));

            n->help(root);

            delete n;
        }

        return root;
    }

    string asString() const {
        HelpNode root;
        return HelpFormatScreen(*asHelpNode(&root));
    }
};


#endif
