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

#include "simulation.h"
#include "postprocessor.h"


REGISTER_SMART_ENUM
(Postprocessor_Factory,
 "Postprocessors are children of Meter objects. The Meter object generates a list of values "
 "and passes them to all the Postprocessors assigned. The Postprocessors themselves perform "
 "operations on that data and pass the new data to their children, which are again Postprocessors. "
 "An Output is a special form of Postprocessor which can not have children but outputs the "
 "data to a file or the screen."
);


//---- Constructors/Destructor ----

Postprocessor::Postprocessor(Node *parent, Simulation *simulation)
	: NodeManyChildren(parent), m_simulation(simulation)
{
}


Postprocessor::~Postprocessor()
{
}



//---- Methods ----

void Postprocessor::setup()
{
  NodeManyChildren::setup();
}


static void do_flush(Node *p) {
    ((Postprocessor*) p)->flush();
}

void Postprocessor::flush()
{
    for_each(m_children.begin(), m_children.end(), do_flush);
}


Node *Postprocessor::instantiateChild(const string &name)
{
    Postprocessor *p = Postprocessor_Factory::byName(name).instantiate(this, m_simulation);
    p->describeInput(&m_output);

    return p;
}

