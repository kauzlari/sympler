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



#include "meter_grouped.h"

//---- Constructors/Destructor ----

MeterGrouped::MeterGrouped(Simulation *simulation)
  : Meter(simulation, 1/*, false*/)
{
    //		cout << "MeterGrouped::MeterGrouped()" << endl;
}


MeterGrouped::MeterGrouped(Simulation *simulation, size_t everyN/*, bool only*/)
  : Meter(simulation, everyN/*, only*/)
{
}


MeterGrouped::~MeterGrouped()
{
}



//---- Methods ----

/*
void MeterGrouped::read(const SGML& e)
{
    Meter::read(e);

    parser p1(e);
    SGML e2;
    while (!p1.eof()) {
        p1.GetSGML(e2);
        if (e2.name == "groups") {
            parser p4(e2);
            set<int>::iterator grIter = groups.begin();
            while (!p4.eof()) {
                grIter = groups.insert(grIter, p4.GetInt());
            }
        }
    }
}
*/


void MeterGrouped::read(const xmlNode *xmln)
{
  Node::read(xmln);

  for (const xmlNode *cur_node = xmln->children; cur_node; cur_node = cur_node->next) {
    if (cur_node->type == XML_ELEMENT_NODE) {
      Node *child;

      child = instantiateChild((const char*) cur_node->name);
      if (child) {
        child->read(cur_node);

        m_children.push_back(child);
      }
    }
  }
//  Node::read(xmln);

/*
  for (const xmlNode *cur_node = xmln; cur_node; cur_node = cur_node->next) {
    if (cur_node->type == XML_ELEMENT_NODE) {
      if (string((const char*) cur_node->name) == "groups") {
        throw gError("MeterGrouped::read", "Groups are not supported right now.");
      } else {
        Node *child;

        child = instantiateChild((const char*) cur_node->name);
        if (child) {
          child->read(cur_node);

          m_children.push_back(child);
        }
      }
    }
  }
*/
}


Node *MeterGrouped::instantiateChild(const string &name)
{
    if (name == "groups")
        return NULL;

    return Meter::instantiateChild(name);
}


