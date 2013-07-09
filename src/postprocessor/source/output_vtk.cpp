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
#include "vtk_support.h"
#include "output_vtk.h"

/* Register this Postprocessor with the factory. */
const Postprocessor_Register<OutputVTK> output_vtk("OutputVTK");


#define VTK_ASCII "ascii"
#define VTK_BINARY "binary"


//---- Constructors/Destructor ----

OutputVTK::OutputVTK(Node *parent, Simulation *simulation)
	: Output(parent, simulation), m_step_counter(0)
{
    init();
}


OutputVTK::~OutputVTK()
{
}



//---- Methods ----

void OutputVTK::init()
{
  m_properties.setClassName("OutputVTK");

  m_properties.setDescription(
    "Outputs the information into a VTK file. This file can be "
    "viewed with ParaView (www.paraview.org)."
  );
    
  STRINGPCOUF
    (name, m_fn,
     "File name of the positions file. This is mandatory.");

  m_properties.addProperty
    ("format", PropertyList::STRING, &m_format,
     new PLCStringList2(VTK_ASCII, VTK_BINARY),
     "Type of the VTK file. Possible values are 'ascii' or 'binary'.");
    
  m_fn = "---";
  m_format = "binary";
}


void OutputVTK::setup()
{
  Output::setup();

//   MSG_DEBUG("OutputVTK::setup","rows=" << m_input_format->rows());

  //    m_idx_positions = m_input_format->indexOf("positions", DataFormat::VECTOR_POINT);
  m_idx_data_format = m_input_format->indexOf(IDDF_DATA_FORMAT, DataFormat::INT);

  if (m_input_format->attrExists(IDDF_N_CELLS))
    m_idx_n_cells = m_input_format->attrByName(IDDF_N_CELLS).index;
  else
    m_idx_n_cells = -1;
}


void OutputVTK::describeInput(DataFormat *input_format)
{
  m_input_format = input_format;
//   MSG_DEBUG("OutputVTK::describeInput","rows=" << m_input_format->rows());

}


void OutputVTK::push(data_sp data)
{
//   MSG_DEBUG("OutputVTK::push","rows=" << m_input_format->rows());


  if (m_fn == "---")
    throw gError("OutputVTK::push", "'nameOutputFile' was not properly defined");

  ofstream s;
  int df = data->intByIndex(m_idx_data_format);
  int_point_t n_cells;
  point_t offset = m_simulation->phase()->boundary()->boundingBox().corner1;
  point_t size = m_simulation->phase()->boundary()->boundingBox().size();
  int idx_start = m_idx_data_format+1;

  s.open(make_filename(m_fn, m_step_counter).c_str());

  s << "# vtk DataFile Version 2.0" << endl;
  s << "SYMPLER simulation " << m_simulation->name() << endl;

  if (m_format == VTK_ASCII)
    s << "ASCII" << endl;
  else
    s << "BINARY" << endl;

  switch (df) {
  case DF_POINTS:
      s << "DATASET UNSTRUCTURED_GRID" << endl;
      break;
  case DF_STRUCTURED_GRID:
      if (m_idx_n_cells < 0 || m_idx_n_cells >= (int) data->format()->rows())
        throw gError
          ("OutputVTK::push",
           "Input data format corrupted. Cannot find number of cells meta information.");

      n_cells = data->intPointByIndex(m_idx_n_cells);

      s << "DATASET STRUCTURED_POINTS" << endl;
      s << "DIMENSIONS " << n_cells.x+1 << " " << n_cells.y+1 << " " << n_cells.z+1 << endl;
      //        s << "ORIGIN 0 0 0 " << endl;
      s << "ORIGIN " << offset.x << " " << offset.y << " " << offset.z << endl;
      s << "SPACING " << size.x/n_cells.x << " " << size.y/n_cells.y << " " << size.z/n_cells.z << endl << endl;
      break;
  default:
      throw gError("OutputVTK", "Don't know how to cast this data into a VTK readable form.");
      break;
  }



  if (df == DF_POINTS) {
    vector<point_t> *p;
    DataFormat::attribute_t attr = data->format()->attrByIndex(idx_start);

    if (attr.datatype != DataFormat::VECTOR_POINT) {
      throw gError
        ("OutputVTK::push", "(Internal error -> check the code) Expecting point "
         "specification for point data format.");
    }

    p = data->vectorPointByIndex(idx_start).value();

//    MSG_DEBUG("OutputVTK::push", "data->format()->rows() = " << data->format()->rows() << ", p->size() = " << p->size());


    s << "POINTS " << p->size() << " double" << endl;

    if (m_format == VTK_ASCII) {
      for (vector<point_t>::iterator i = p->begin(); i != p->end(); i++) {
        s << i->x << " " << i->y << " " << i->z << endl;
      }
    } else {
      for (vector<point_t>::iterator i = p->begin(); i != p->end(); i++) {
        point_t p = *i;
        vtkSwapPoint(p);
        s.write((char*) &p.x, sizeof(double));
        s.write((char*) &p.y, sizeof(double));
        s.write((char*) &p.z, sizeof(double));
      }
    }

    s << endl;

    s << "CELLS " << p->size() << " " << 2*p->size() << endl;

    if (m_format == VTK_ASCII) {
      for (size_t i = 0; i < p->size(); ++i) {
        s << "1 " << i << endl;
      }
    } else {
      for (size_t i = 0; i < p->size(); ++i) {
        int h = 1, j = i;
        vtkSwapInt(h);
        vtkSwapInt(j);
        s.write((char*) &h, sizeof(int));
        s.write((char*) &j, sizeof(int));
      }
    }

    s << endl;

    s << "CELL_TYPES " << p->size() << endl;

    if (m_format == VTK_ASCII) {
      for (size_t i = 0; i < p->size(); ++i) {
        s << "1" << endl;
      }
    } else {
      for (size_t i = 0; i < p->size(); ++i) {
        int h = 1;
        vtkSwapInt(h);
        s.write((char*) &h, sizeof(int));
      }            
    }

    s << endl;

    s << "POINT_DATA " << p->size() << endl;

    idx_start++;
  } else
    s << "CELL_DATA " << n_cells.x*n_cells.y*n_cells.z << endl;

//  MSG_DEBUG("OutputVTK::push", "idx_start = " << idx_start << ", rows = " << m_input_format->rows());

  for (size_t i = idx_start; i < m_input_format->rows(); ++i) {
    DataFormat::attribute_t attr = m_input_format->attrByIndex(i);

//             MSG_DEBUG("OutputVTK::push", "attr.name = " << attr.name << ", attr.index = " << attr.index);

    switch (attr.datatype) {
    case DataFormat::VECTOR_INT:
        if (m_format == VTK_ASCII)
          vtkWriteVectorInt_Ascii(s, attr.name, data->vectorIntByIndex(i));
        else
          vtkWriteVectorInt_Binary(s, attr.name, data->vectorIntByIndex(i));
        break;
    case DataFormat::VECTOR_DOUBLE:
        if (m_format == VTK_ASCII)
          vtkWriteVectorDouble_Ascii(s, attr.name, data->vectorDoubleByIndex(i));
        else
          vtkWriteVectorDouble_Binary(s, attr.name, data->vectorDoubleByIndex(i));
        break;
    case DataFormat::VECTOR_POINT:
        if (m_format == VTK_ASCII)
          vtkWriteVectorPoint_Ascii(s, attr.name, data->vectorPointByIndex(i));
        else
          vtkWriteVectorPoint_Binary(s, attr.name, data->vectorPointByIndex(i));
        break;
    case DataFormat::VECTOR_TENSOR:
        if (m_format == VTK_ASCII)
          vtkWriteVectorTensor_Ascii(s, attr.name, data->vectorTensorByIndex(i));
        else
          vtkWriteVectorTensor_Binary(s, attr.name, data->vectorTensorByIndex(i));
        break;
    default:
        throw gError
          ("OutputVTK::push",
           "Unsupported data format for attribute '" + attr.name + "'");
    }

    s << endl;
  }

  s.close();

  m_step_counter++;

}



