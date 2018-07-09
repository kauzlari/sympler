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
#include "output_mathematica.h"

/* Register this Postprocessor with the factory. */
const Postprocessor_Register<OutputMathematica> output_mathematica("OutputMathematica");



//---- Constructors/Destructor ----

OutputMathematica::OutputMathematica(Node *parent, Simulation *simulation)
  : OutputFile(parent, simulation), m_comma(false)
{
    init();
}


OutputMathematica::~OutputMathematica()
{
}

void OutputMathematica::init()
{
  m_properties.setClassName("OutputMathematica");

  m_properties.setDescription(
    "Output the information obtained from the preceding postprocessor "
    "queue into a Mathematica readable file."
  );

  m_comment_start = "(*";
  m_comment_end = "*)";
}


//---- Methods ----

void OutputMathematica::push(data_sp data)
{
    openStream();
  
    if (!m_multiple_files) {
        if (m_comma)
            m_s << ", ";
        else
            m_s << "{";

        m_comma = true;
    }

    m_s << "{";

    bool comma = false;
    for (list<int>::iterator i = m_columns.begin(); i != m_columns.end(); i++) {
        if (comma)
            m_s << ", ";
        comma = true;
        m_s << data->toMathematicaByIndex(*i, (m_multiple_files?"\n":""));
    }
  
    m_s << "}";

    closeStream();
}


void OutputMathematica::flush()
{
  if(!m_multiple_files) {
      m_s << "}";
      
      OutputFile::flush();
  }
}

