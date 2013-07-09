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
#include "output_pdb.h"

/* Register this Postprocessor with the factory. */
const Postprocessor_Register<OutputPDB> output_pos_vel_pdb("OutputPDB");


//---- Constructors/Destructor ----

OutputPDB::OutputPDB(Node *parent, Simulation *simulation)
	: Output(parent, simulation), m_step_counter(0)
{
    init();
}


OutputPDB::~OutputPDB()
{
}



//---- Methods ----

void OutputPDB::init()
{
  m_properties.setClassName("OutputPDB");

  m_properties.setDescription
    ("Writes position information into a file using the PDB format."
     );
    
  STRINGPCOUF
    (namePositions, m_posfn,
     "File name of the positions file. This is mandatory.");
  m_posfn = "---";

  BOOLPC
    (multipleFiles, m_multiple_files,
     "Create a seperate file for each timestep.");
    
  m_posfn = "positions.out";
  m_multiple_files = false;
}


void OutputPDB::setup()
{
  Output::setup();
	
  if(m_posfn == "---") 
    throw gError("OutputPDB::read: 'positionsfn' not properly defined.");
  if (!m_multiple_files) {
    m_pos_s.open(m_posfn.c_str());
  }
}


void OutputPDB::describeInput(DataFormat *input_format)
{
    m_input_format = input_format;
	
    m_idx_positions = m_input_format->indexOf("positions", DataFormat::VECTOR_POINT);
}



string OutputPDB::point2string(const point_t &p, bool comma)
{
	stringstream s;
	s.flags(ios::fixed);

	for (int i = 0; i < SPACE_DIMS; i++) {
		if (i) {
			if (comma)
				s << ", ";
			else
				s << " ";
		}
		s << setiosflags(ios::fixed) << setprecision(3) << setw(7) << p[i];
	}
	
	return s.str();
}

void OutputPDB::writeVectorPoint(ostream &s, vector_point_sp vp)
{
    int line = 0;

	if (m_step_counter && !m_multiple_files)
		s << endl;
		
	s << "COMPND    particle positions" << endl;

    for (vector<point_t>::iterator i = vp->begin(); i != vp->end(); i++) {
        s << "HETATM " << setw(4) << (line+1) << "  H           1     " << point2string(*i, false) << "  1.00  0.00           H" << endl;
		
        line++;
    }
	
	s << "TER " << setw(7) << (line+1) << endl;
	s << "END" << endl;
}


void OutputPDB::push(data_sp data)
{
	if (m_multiple_files) {
		m_pos_s.open(make_filename(m_posfn, m_step_counter).c_str());
	}

	writeVectorPoint(m_pos_s, data->vectorPointByIndex(m_idx_positions));
    //	writeVectorPoint(m_vel_s, data->vectorPointByIndex(m_idx_velocities));
	++m_step_counter;
	
	if (m_multiple_files) {
		m_pos_s.close();
        //		m_vel_s.close();
	}
}


void OutputPDB::flush()
{
	if (!m_multiple_files) {
		m_pos_s.close();
        //		m_vel_s.close();
	}
}

