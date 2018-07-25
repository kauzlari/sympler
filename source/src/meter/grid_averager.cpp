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



#include "grid_averager.h"

#include "grid_meter.h"
#include "simulation.h"

using namespace std;


//---- Constructors/Destructor ----

GridAverager::GridAverager(Simulation *simulation)
  : Meter(simulation, 1), m_step(0)
{
    init();
}


GridAverager::~GridAverager()
{
}



//---- Methods ----

void GridAverager::init()
{
  m_properties.setClassName("GridAverager");
    
  INTPC
    (avgOver, m_avg_steps, 0,
     "Number of time steps to average over. Note that this may be at maximum equal to 'measureEvery'. The difference 'measureEvery'-'avgOver' represents the time interval, in which this GridAverager does nothing.");

  m_avg_steps = 1;
}


void GridAverager::setup()
{
  if(m_avg_steps > m_measure_every_n) 
    throw gError("GridAverager::setup()", "'avgOver' > 'measureEvery' does not make any sense");
  
  size_t nOfC = ((Simulation*) m_parent)->phase()->nColours();
  // should create empty vector<int> in every slot
  m_locations.resize(nOfC);   
  m_n_particles.resize(nOfC);


//   MSG_DEBUG("GridAverager::setup", "setting up my GridMeters");
  for(list<GridMeter*>::iterator gm = m_meters.begin(); gm != m_meters.end(); ++gm)
    {
      MSG_DEBUG("GridAverager::setup", "Starting setup for GridMeter \"" << (*gm)->className() << "\".");
      (*gm)->setup();
    }
  /* Make sure, setup of the postprocessor is called *last*. */
  Meter::setup();
}


Node *GridAverager::instantiateChildWithXML(const string &name, const xmlNode *xmln)
{
  if (GridMeter_Factory::exists(name)) {
    GridMeter *m = GridMeter_Factory::byName(name).instantiate(this);

    m->read(xmln);

    m_meters.push_back(m);
  } else
    return Meter::instantiateChild(name);

  return NULL;
}


void GridAverager::average(data_sp data)
{
    for (list<GridMeter*>::iterator i = m_meters.begin(); i != m_meters.end(); i++) {
        (*i)->measure(data);
    }
}


void GridAverager::finishStep(data_sp data)
{
//   MSG_DEBUG("GridAverager::finishStep", "finishing: m_step = " << m_step);

  for (list<GridMeter*>::iterator i = m_meters.begin(); i != m_meters.end(); i++) {
    (*i)->finishStep(data, m_step);
  }

  //  MSG_DEBUG("GridAverager::finishStep", "averaging");

  // below, we divide through m_step so that also the last INCOMPLETE average at the end
  // of the simulation is computed correctly

  if (m_avg_steps != 1) {
    for (size_t i = m_idx_data_start; i < m_format.rows(); i++) {
      DataFormat::attribute_t attr;

      attr = m_format.attrByIndex(i);

      //      MSG_DEBUG("GridAverager::finishStep", "attr.name = " << attr.name);
      //      MSG_DEBUG("GridAverager::finishStep", "m_step = " << m_step);

      switch (attr.datatype) {
      case DataFormat::VECTOR_DOUBLE: {
        vector<double> *v = data->vectorDoubleByIndex(i).value();
        for (vector<double>::iterator j = v->begin(); j != v->end(); j++)
          *j /= /*m_avg_steps*/ m_step;
        break;
      }
      case DataFormat::VECTOR_POINT: {
        vector<point_t> *v = data->vectorPointByIndex(i).value();
        for (vector<point_t>::iterator j = v->begin(); j != v->end(); j++)
          *j /= /*m_avg_steps*/ m_step;
        break;
      }
      default:
          throw gError
            ("GridAverager::finishStep", "Unsupported data "
             "format for attribute '" + attr.name + "'");
      }
    }
  }

//  MSG_DEBUG("GridAverager::finishStep", "distributing");

  distribute(data);
}


void GridAverager::writeMoreBody(ostream &s, int shift)
{
    for (list<GridMeter*>::iterator i = m_meters.begin(); i != m_meters.end(); i++)
        (*i)->write(s, shift+1);
}


