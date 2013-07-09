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



#include <stdlib.h>

#include <string>
#include <ostream>
#include <iomanip>

#include "simulation.h"
#include "output_dcd.h"
#include "vtk_support.h"

/* Register this Postprocessor with the factory. */
const Postprocessor_Register<OutputDCD> output_dcd("OutputDCD");



//---- Constructors/Destructor ----

OutputDCD::OutputDCD(Node *parent, Simulation *simulation)
	: Output(parent, simulation), m_step_counter(0)
{
    init();
}


OutputDCD::~OutputDCD()
{
}



//---- Methods ----

void OutputDCD::init()
{
  m_properties.setClassName("OutputDCD");

  m_properties.setDescription(
    "Outputs the information into a DCD file. This file can be "
    "viewed with VMD (http://www.ks.uiuc.edu/Research/vmd/). Note "
    "that this data format CANNOT store additional information, such "
    "as, e.g., internal energy. Use OutputVTK instead."
  );
    
  STRINGPCOUF
    (name, m_fn,
     "File name of the positions file.");

  m_fn = "trajectory.dcd";
}


void OutputDCD::setup()
{
  Output::setup();
  m_idx_data_format = m_input_format->indexOf(IDDF_DATA_FORMAT, DataFormat::INT);
}


void OutputDCD::describeInput(DataFormat *input_format)
{
  m_input_format = input_format;
}


void OutputDCD::push(data_sp data)
{
  if (m_fn == "---")
    throw gError("OutputDCD::push", "'nameOutputFile' was not properly defined");

  int df = data->intByIndex(m_idx_data_format);
  point_t offset = m_simulation->phase()->boundary()->boundingBox().corner1;
  point_t size = m_simulation->phase()->boundary()->boundingBox().size();
  int idx_start = m_idx_data_format+1;

  if (df != DF_POINTS) {
    throw gError
      ("OuputDCD::push",
       "The DCD format can only hold trajectory data. Please do only use MeterPosVel "
       "with this output module.");
  }

  vector<point_t> *p;
  DataFormat::attribute_t attr = data->format()->attrByIndex(idx_start);

  if (attr.datatype != DataFormat::VECTOR_POINT) {
    throw gError
      ("OutputDCD::push", "(Internal error) Expecting point "
       "specification for point data format.");
  }

  p = data->vectorPointByIndex(idx_start).value();

  if (!m_s.is_open()) {
    int h;
    double d;
//     float f;
    char buf[100];
    char head[] = "CORD";
    char title[80] =
      "SYMPLER trajectory file                                                        ";

    m_s.open(m_fn.c_str());

    h = 84; /* Magic number ????? */
    m_s.write((char*) &h, sizeof(int));

    /* Header */
    m_s.write((char*) head, 4*sizeof(char));

    h = m_simulation->controller()->timesteps(); /* Number of frames */
    m_s.write((char*) &h, sizeof(int));

    h = 0; /* Number of steps in previous runs */
    m_s.write((char*) &h, sizeof(int));

    h = 1; /* Step interval */
    m_s.write((char*) &h, sizeof(int));

    h = 0; /* Number of total steps */
    m_s.write((char*) &h, sizeof(int));

    memset(buf, 0, 3*sizeof(int)); /* 3 zeros */
    m_s.write((char*) buf, 3*sizeof(int));

    h = 0; /* Number of degrees of freedom, zero for now??? */
    m_s.write((char*) &h, sizeof(int));

    h = 0; /* Number of frozen atoms */
    m_s.write((char*) &h, sizeof(int));

    d = m_simulation->controller()->dt(); /* Time step */
    m_s.write((char*) &d, sizeof(double));

    //    h = 0; /* Crystallographic group */
    //    m_s.write((char*) &h, sizeof(int));

    memset(buf, 0, 8*sizeof(int)); /* 8 zeros */
    m_s.write((char*) buf, 8*sizeof(int));

    h = 0; /* Version number */
    m_s.write((char*) &h, sizeof(int));

    h = 84; /* Magic number ????? */
    m_s.write((char*) &h, sizeof(int));

    h = 84; /* Title size */
    m_s.write((char*) &h, sizeof(int));

    h = 1; /* Number of lines */
    m_s.write((char*) &h, sizeof(int));

    m_s.write((char*) title, 80*sizeof(char));

    h = 84; /* Title size, again */
    m_s.write((char*) &h, sizeof(int));

    h = 4; /* No idea what this could be */
    m_s.write((char*) &h, sizeof(int));

    m_n_atoms = p->size(); /* Number of atoms */
    m_s.write((char*) &m_n_atoms, sizeof(int));

    h = 4; /* No idea what this could be */
    m_s.write((char*) &h, sizeof(int));

    /* That's it! */
  }


  if (p->size() != m_n_atoms)
    throw gError
      ("OutputDCD::push",
       "The number of particles is not allowed to change when using the DCD file format.");

  int out_integer = m_n_atoms*sizeof(int);


  m_s.write((char*) &out_integer, sizeof(int));
  for (vector<point_t>::iterator i = p->begin(); i != p->end(); i++) {
    point_t p = *i;
    float f;
    //    vtkSwapPoint(p);
    f = p.x;
    m_s.write((char*) &f, sizeof(float));
  }
  m_s.write((char*) &out_integer, sizeof(int));

  m_s.write((char*) &out_integer, sizeof(int));
  for (vector<point_t>::iterator i = p->begin(); i != p->end(); i++) {
    point_t p = *i;
    float f;
    //    vtkSwapPoint(p);
    f = p.y;
    m_s.write((char*) &f, sizeof(float));
  }
  m_s.write((char*) &out_integer, sizeof(int));

  m_s.write((char*) &out_integer, sizeof(int));
  for (vector<point_t>::iterator i = p->begin(); i != p->end(); i++) {
    point_t p = *i;
    float f;
    //    vtkSwapPoint(p);
    f = p.z;
    m_s.write((char*) &f, sizeof(float));
  }
  m_s.write((char*) &out_integer, sizeof(int));

  m_step_counter++;
}


void OutputDCD::flush()
{    
  m_s.close();
}

