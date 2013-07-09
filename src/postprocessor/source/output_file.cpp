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

#include <ostream>
#include <algorithm>

#include "simulation.h"
#include "output_file.h"


#define M_METER ((Meter*) m_parent)
#define M_SIMULATION  ((Simulation*)(M_METER->parent()))
#define M_CONTROLLER M_SIMULATION->controller()


/* Register this Postprocessor with the factory. */
const Postprocessor_Register<OutputFile> output_file("OutputFile");


//---- Constructors/Destructor ----

OutputFile::OutputFile(Node *parent, Simulation *simulation)
	: Output(parent, simulation), m_step_counter(0)
{
    init();
}


OutputFile::~OutputFile()
{
}



//---- Methods ----


void OutputFile::init()
{
  m_properties.setClassName("OutputFile");

  m_properties.setDescription(
    "Output the information obtained from the preceding postprocessor "
    "queue into a text file separated by spaces.\n"
    "\n"
    "The attribute 'multipleFiles' can be used to output a file for each timestep. If multipleFiles = \"yes\", You are allowed to output either vector-data or single-valued data but no mix of both. Note that, If multipleFiles = \"no\", output of array data is usually not practical since all data will be stored in a single line per output-event. Array-data is for example non-averaged data for each particle or data from a GridAverager. MeterAverage is a GridAverager too. But since it produces array-data of length 1, here, multipleFiles = \"yes\" makes sense."
  );
    
  STRINGPCOUF
    (name, m_filename,
     "Output filename.");
  STRINGPC
    (columns, m_cols,
     "Columns to write to disk. The different columns have to be separated by the pipe "
     "symbol_'|'. The information can "
     "be inverted. Using an exclamation mark as the first character outputs all columns "
     "*except* for the specified in this attribute (Example: '!time').");
    
  BOOLPC
    (writeHead, m_do_head,
     "Write a header stating the column names to the first line of the file.");

  BOOLPC
    (multipleFiles, m_multiple_files,
     "Create a separate file for each timestep.");
    
  BOOLPC
    (binary, m_binary,
     "If set to \"no\", ascii files are produced, if set to \"yes\", binary files. In binary mode, of course no header is written.");
    
  m_filename = m_cols = "---";
  m_do_head = true;
  m_multiple_files = true;
  m_comment_start = "#";
  m_comment_end = "";
  m_binary = false;
}


void OutputFile::setup()
{
  Output::setup();

  M_CONTROLLER->registerForSetupAfterParticleCreation(this);
}

void OutputFile::setupAfterParticleCreation()
{

  if (m_filename == "---") {
    throw gError
      ("OutputFile::read",
       "No file name given.");
  }

  if (m_cols == "---") {
    MSG_DEBUG("OutputFile::setupAfterParticleCreation", "for " << m_parent->className() << ": user selection for columns: ---");
    for (size_t i = 0; i < m_input_format->rows(); ++i)
      m_columns.push_back(i);
  } 
  else {
    MSG_DEBUG("OutputFile::setupAfterParticleCreation", "for " << m_parent->className() << ": user selection for columns: " << m_cols);
    string cols = m_cols;
    bool run = true;
    bool invert = false;
        
    /* If the first character is an exclamation mark, all the columns
       given in 'columns' are EXCLUDED.
    */
    if (cols[0] == '!') {
      invert = true;

      cols = string(m_cols, 1);

      for (size_t i = 0; i < m_input_format->rows(); ++i) {
	       m_columns.push_back(i);
      }
    }

    while (run) {
      string cur;
      size_t pos = cols.find('|');

      if (pos == string::npos) {
	       run = false;
	       cur = cols;
      } 
      else {
	       cur = string(cols, 0, pos);
	       cols = string(cols, pos+1);
      }

      if (m_input_format->attrExists(cur)) {
	if (invert)
	  m_columns.remove(m_input_format->attrByName(cur).index);
	else
	  m_columns.push_back(m_input_format->attrByName(cur).index);
      } 
      else
	throw gError
	  ("OutputFile::setup",
           "For " + m_parent->className() + ": No column '" + cur + "' in output stream. "
	   "Possibilities are: " + m_input_format->toString());
    }
  } // end of else of if (m_cols == "---")

  /* If we use multiple files we can as well output vector_doubles correctly. */
  m_all_vectors = false;

  if (m_multiple_files) {
    m_all_vectors = true;
    bool one_vec = false;
    bool one_single = false;
    
    for (list<int>::iterator i = m_columns.begin(); i != m_columns.end(); i++) {
      DataFormat::attribute_t attr = m_input_format->attrByIndex(*i);

      if (attr.datatype != DataFormat::VECTOR_DOUBLE && attr.datatype != DataFormat::VECTOR_INT && attr.datatype != DataFormat::VECTOR_POINT)
      {
        MSG_DEBUG("OutputFile::setup",
                  "For " << m_parent->className() << ": multiple files + single: " << m_input_format->attrByIndex(*i).name);
        m_all_vectors = false;
        one_single = true;
        if(one_vec) // ...then abort due to mixing of data
        {
          string single_data, array_data;
          for (list<int>::iterator i_err = m_columns.begin(); i_err != m_columns.end(); i_err++) {
            DataFormat::attribute_t attr_err = m_input_format->attrByIndex(*i_err);

            if (attr_err.datatype != DataFormat::VECTOR_DOUBLE && attr_err.datatype != DataFormat::VECTOR_INT && attr_err.datatype != DataFormat::VECTOR_POINT)
            {
              if(single_data != "")
                single_data += ", ";
              single_data += m_input_format->attrByIndex(*i_err).name;
            }
            else
            {
              if(array_data != "")
                array_data += ", ";
              array_data += m_input_format->attrByIndex(*i_err).name;
            }
          }
          throw gError("OutputFile::setup",
                       "For " + m_parent->className() + ": For multipleFiles = \"yes\", don't mix array-data and single data in one output file. Probably you should check the settings for attribute \"columns\".\nArray-data you told me to write: " + array_data + "\nSingle data you told me to write: " + single_data);
        }
        
      } // end: if(no vector)
      else
      {
        MSG_DEBUG("OutputFile::setup",
                  "For " << m_parent->className() << ": multiple files + array: " << m_input_format->attrByIndex(*i).name);
        one_vec = true; 
        if(one_single) // ...then abort due to mixing of data
        {
          string single_data;
          string array_data;
          for (list<int>::iterator i_err = m_columns.begin(); i_err != m_columns.end(); i_err++) {
            DataFormat::attribute_t attr_err = m_input_format->attrByIndex(*i_err);
  
            if (attr_err.datatype != DataFormat::VECTOR_DOUBLE && attr_err.datatype != DataFormat::VECTOR_INT && attr_err.datatype != DataFormat::VECTOR_POINT)
            {
              MSG_DEBUG("OutputFile::setup", "single_data = " << single_data);
              if(single_data != "")
                single_data += ", ";
              single_data += m_input_format->attrByIndex(*i_err).name;
              MSG_DEBUG("OutputFile::setup", "single_data = " << single_data);
            }
            else
            {
              MSG_DEBUG("OutputFile::setup", "array_data = " << array_data);
              if(array_data != "")
                array_data += ", ";
              array_data += m_input_format->attrByIndex(*i_err).name;
            }
          }
          throw gError("OutputFile::setup",
                       "For " + m_parent->className() + ": For multipleFiles = \"yes\", don't mix array-data and single data in one output file. Probably you should check the settings for attribute \"columns\".\nArray-data you told me to write: " + array_data + "\nSingle data you told me to write: " + single_data);
        }
      } // end: else of if(no vector)
    } // end: for-loop over m_columns
  } // end: if(m_multiple_files)
  else {
    for (list<int>::iterator i = m_columns.begin(); i != m_columns.end(); i++) {
      DataFormat::attribute_t attr = m_input_format->attrByIndex(*i);

      if (attr.datatype == DataFormat::VECTOR_DOUBLE || attr.datatype == DataFormat::VECTOR_INT || attr.datatype == DataFormat::VECTOR_POINT)
      {
        string array_data;
        // report error because no vectors allowed in single file
        for (list<int>::iterator i_err = m_columns.begin(); i_err != m_columns.end(); i_err++) {
          DataFormat::attribute_t attr_err = m_input_format->attrByIndex(*i_err);
  
          if (attr_err.datatype == DataFormat::VECTOR_DOUBLE || attr_err.datatype == DataFormat::VECTOR_INT || attr_err.datatype == DataFormat::VECTOR_POINT)
          {
            if(array_data != "")
              array_data += ", ";
            array_data += m_input_format->attrByIndex(*i_err).name;
          }
        }
        // currently commented out because I want to write those with length one (e.g. from MeterAverage)        
/*        throw gError("OutputFile::setup",
                     "For " + m_parent->className() + ": For multipleFiles = \"no\", I cannot write array-data.\nArray-data you told me to write: " + array_data);*/
      }
    } // end: loop over m_columns   
    // not boiled out? so open the single file
    MSG_DEBUG("OutputFile::setup", "Opening single file " << m_filename.c_str());
    if(m_binary) {
      m_s.open(m_filename.c_str(), ios::out|ios::binary); 
    }     
    else {
      m_s.open(m_filename.c_str());
      m_s.flags(ios::scientific);
      writeHeader();
    }

  } // end: else of if(m_multiple_files)

}


void OutputFile::describeInput(DataFormat *input_format)
{
  /* Nothing to complain about, we're just outputting everything. */
  m_input_format = input_format;

}


void OutputFile::push(data_sp data)
{
  openStream();

  if(m_binary) {
    size_t size;
    DataFormat::attribute_t attr = m_input_format->attrByIndex(*m_columns.begin());
    
    if (attr.datatype == DataFormat::VECTOR_DOUBLE)
      size = data->vectorDoubleByIndex(*m_columns.begin())->size();
    else if (attr.datatype == DataFormat::VECTOR_INT)
      size = data->vectorIntByIndex(*m_columns.begin())->size();
    else if (attr.datatype == DataFormat::VECTOR_POINT)
      size = data->vectorPointByIndex(*m_columns.begin())->size();
    else
      throw gError("OutputFile::push", "Unsupported datatype for binary output for attribute \"" + attr.name + "\".");
    
    for (size_t i = 0; i < size; i++) {
      for (list<int>::iterator j = m_columns.begin(); j != m_columns.end(); j++) {
	attr = m_input_format->attrByIndex(*j);
	
	if (attr.datatype == DataFormat::VECTOR_DOUBLE) {
	  double d = (*data->vectorDoubleByIndex(*j))[i];
//   	  MSG_DEBUG("OutputFile::push", "VECTOR_DOUBLE case: i=" << i << ", size=" << size << ", d=" << d);
	  // vtkSwapDouble(d);
	  m_s.write((char*) &d, sizeof(double));

//  	  MSG_DEBUG("OutputFile::push", "flushing now! (remove for real sim!)");
//  	  m_s.flush();
	}
	else if (attr.datatype == DataFormat::VECTOR_INT) {
	  int num = (*data->vectorIntByIndex(*j))[i];
	  // vtkSwapInt(num);
	  m_s.write((char*) &num, sizeof(int));
	}
	else if (attr.datatype == DataFormat::VECTOR_POINT) {

// 	  MSG_DEBUG("OutputFile::push", "VECTOR_POINT case: i=" << i);

	  point_t p = (*data->vectorPointByIndex(*j))[i];
	  // 	    double d = p.x;
	  // vtkSwapDouble(d);

	  //---START: DEBUGGING--------------------
// 	  ifstream debugStream;
// 	  ifstream::pos_type dbgsize;

// 	  debugStream.open("interpoldat/dataPN_NS_00001.bin", ios::in|ios::binary|ios::ate);
// 	  dbgsize = debugStream.tellg();
// 	  MSG_DEBUG("OutputFile::push", "size before: " << dbgsize);
// 	  debugStream.close();

// 	  MSG_DEBUG("OutputFile::push", "writing " << p.x);
	  //---END: DEBUGGING--------------------

	  m_s.write((char*) &(p.x), sizeof(double));

	  //---START: DEBUGGING--------------------
// 	  debugStream.open("interpoldat/dataPN_NS_00001.bin", ios::in|ios::binary|ios::ate);
// 	  if(debugStream.is_open()) {
// 	    dbgsize = debugStream.tellg();
// 	    MSG_DEBUG("OutputFile::push", "size after: " << dbgsize);
// 	    debugStream.close();
// 	  }
// 	  else MSG_DEBUG("OutputFile::push", "UNABLE TO OPEN");
	  //---END: DEBUGGING--------------------

	  // 	    double d = p.y;
	  // vtkSwapDouble(d);
	  m_s.write((char*) &(p.y), sizeof(double));


	  //---START: DEBUGGING--------------------
// 	  MSG_DEBUG("OutputFile::push", "just written now flushing" << p.y);
// 	  m_s.flush();

// 	  debugStream.open("interpoldat/dataPN_NS_00001.bin", ios::in|ios::binary|ios::ate);
// 	  dbgsize = debugStream.tellg();
// 	  MSG_DEBUG("OutputFile::push", "size after: " << dbgsize);
// 	  debugStream.close();
	  //---END: DEBUGGING--------------------

	  // 	    double d = p.z;
	  // vtkSwapDouble(d);
	  m_s.write((char*) &(p.z), sizeof(double));

	  //---START: DEBUGGING--------------------
// 	  MSG_DEBUG("OutputFile::push", "just written now flushing" << p.z);
// 	  m_s.flush();

// 	  debugStream.open("interpoldat/dataPN_NS_00001.bin", ios::in|ios::binary|ios::ate);
// 	  dbgsize = debugStream.tellg();
// 	  MSG_DEBUG("OutputFile::push", "size after: " << dbgsize);
// 	  debugStream.close();
	  //---END: DEBUGGING--------------------

	}
      } // (for(list<int>::iterator ...)
    } // end of for(size_t i ...
  } // end of if(m_binary)
  else {
    if (m_all_vectors) {
      size_t size;
      DataFormat::attribute_t attr = m_input_format->attrByIndex(*m_columns.begin());
      
      if (attr.datatype == DataFormat::VECTOR_DOUBLE)
	size = data->vectorDoubleByIndex(*m_columns.begin())->size();
      else if (attr.datatype == DataFormat::VECTOR_INT)
	size = data->vectorIntByIndex(*m_columns.begin())->size();
      else if (attr.datatype == DataFormat::VECTOR_POINT)
	size = data->vectorPointByIndex(*m_columns.begin())->size();
      else
	throw gError("OutputFile::push", "Unsupported datatype.");
      
      for (size_t i = 0; i < size; i++) {
	for (list<int>::iterator j = m_columns.begin(); j != m_columns.end(); j++) {
	  attr = m_input_format->attrByIndex(*j);
	  
	  if (attr.datatype == DataFormat::VECTOR_DOUBLE) {
	    m_s << (*data->vectorDoubleByIndex(*j))[i] << " ";
	  }
	  else if (attr.datatype == DataFormat::VECTOR_INT) {
	    m_s << (*data->vectorIntByIndex(*j))[i] << " ";
	  }
	  else if (attr.datatype == DataFormat::VECTOR_POINT) {
	    point_t p = (*data->vectorPointByIndex(*j))[i];
	    
	    m_s << p.x << " " << p.y << " " << p.z << " ";
	  }
	}
	m_s << endl;
      }
    } // end if(m_all_vectors) 
    else {
      for (list<int>::iterator i = m_columns.begin(); i != m_columns.end(); i++) {
	
	//       MSG_DEBUG("OutputFile::push", "Now format flags = " << m_s.flags());
	//       MSG_DEBUG("OutputFile::push", (*i));
	m_s << data->toStringByIndex(*i, (m_multiple_files?"\n":"")) << " ";
      }
      m_s << endl; 
    }
  } // end of else of if if(m_binary)
  closeStream();
}


void OutputFile::flush()
{    
    if (!m_multiple_files)
        m_s.close();
}


//---


void OutputFile::writeHeader()
{
  if(m_do_head) {
    m_s << m_comment_start;
    for (list<int>::iterator i = m_columns.begin(); i != m_columns.end(); i++)
      m_s << m_input_format->attrByIndex(*i).name << " ";
    m_s << m_comment_end << endl;
  }
}


void OutputFile::openStream()
{
  if(m_multiple_files) {    
    if(m_binary) {
//       MSG_DEBUG("OutputFile::openStream", "multiple-binary case");
      m_s.open(make_filename(m_filename, m_step_counter).c_str(), ios::out|ios::binary);      
    }
    else {
      m_s.open(make_filename(m_filename, m_step_counter).c_str());
      writeHeader();
    }
    ++m_step_counter;      
  }

}


void OutputFile::closeStream()
{
    if (m_multiple_files) {
        m_s.close();
    }
}

