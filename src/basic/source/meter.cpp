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



#include <string>
#include <algorithm>

#include "simulation.h"
#include "postprocessor.h"

using namespace std;

REGISTER_SMART_ENUM
(Meter_Factory,
 "Meters are object that determine quantities from a running simulation "
 "and pass them to a postprocessor object (see Postprocessors).");


#define M_SIMULATION ((Simulation*) m_parent)


//---- Constructors/Destructor ----

Meter::Meter(Simulation *simulation)
    : NodeManyChildren((Node*) simulation), 
      m_measure_every_n(1), m_counter(1)
{
    init();
}


Meter::Meter(Simulation *simulation, size_t everyN)
    : NodeManyChildren((Node*) simulation), m_measure_every_n(everyN),
      m_counter(everyN)
{
    init();
}


/*virtual*/ Meter::~Meter()
{
}


void Meter::init()
{
  m_properties.setClassName("Meter");
    
  INTPC
    (measureEvery, m_measure_every_n, 0,
     "Number of timesteps between measurements.");
  INTPC
    (fromStepOn, m_from_step_on, -1,
     "Step at which measurements are started.");
    
  m_measure_every_n = 1;
  m_from_step_on = 0;
}


void Meter::aboutToStart()
{
}


/*virtual*/ void Meter::setup()
{
  NodeManyChildren::setup();

  m_counter = m_measure_every_n;
  m_from_time_on = (m_from_step_on)*M_SIMULATION->controller()->dt();
}

Node *Meter::instantiateChild(const string &name)
{
    Postprocessor *p = Postprocessor_Factory::byName(name).instantiate(this, (Simulation*) m_parent);
    p->describeInput(&m_format);
//     MSG_DEBUG("Meter::instantiateChild", "m_format.toString() = " << m_format.toString());
    return p;
}


static void do_flush(Node *p) {
    ((Postprocessor*) p)->flush();
}

void Meter::flush()
{
    for_each(m_children.begin(), m_children.end(), do_flush);
}
