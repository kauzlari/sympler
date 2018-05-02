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


#include "apply_vector_field_file.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


const Callable_Register<ApplyVectorFieldFile> apply_vector_field_file("ApplyVectorFieldFile");

/* --- ApplyVectorFieldFile --- */

ApplyVectorFieldFile::ApplyVectorFieldFile(Simulation* sim)
  : Thermostat(sim), m_colour(11111111)
{
  init();
}


void ApplyVectorFieldFile::init()
{
  m_properties.setClassName("ApplyVectorFieldFile");

  m_properties.setDescription(
    "When called, this callable adds a user-specified vector-field (from a binary file) to each particle. In terms of rows (slow index) and columns (fast index) the binary file must contain the cartesian components as columns and the particles as rows. "
  );

  STRINGPC
    (species, m_species,
     "Species this callable should work on.");

  m_species = "UNDEF";

  STRINGPC
    (symbol, m_symbolName,
     "Name of the symbol to be modified.");

  m_symbolName = "undefined";

  STRINGPCINF
    (name, m_fileName,
     "Name of the binary file containing the vector field.");

  m_fileName = "undefined";

}


void ApplyVectorFieldFile::setup()
{
  Thermostat::setup();

  if(m_species == "UNDEF")
    throw gError("ApplyVectorFieldFile::setup", "Attribute 'species' was not defined!");

  m_colour = M_MANAGER->getColour(m_species);

  if(m_symbolName == "undefined")
    throw gError("ApplyVectorFieldFile::setup", "Attribute 'symbol' was not defined!");

  // the attribute should already exist; the following should also search and check the format
  try
    {
      m_offset = Particle::s_tag_format[m_colour].indexOf(m_symbolName, DataFormat::POINT);
      m_offset = Particle::s_tag_format[m_colour].offsetByIndex(m_offset);
    }
  catch(gError& err)
    {
      throw gError("ApplyVectorFieldFile::setup", "search for symbol failed. The message was " + err.message()); 
    }
  
  if(m_fileName == "undefined")
    throw gError("ApplyVectorFieldFile::setup", "Attribute 'fileName' was not defined!");

  ifstream inFile;
  ifstream::pos_type fileSize;
  inFile.open(m_fileName.c_str(), ios::in|ios::binary|ios::ate);

  if(inFile.is_open()) {
    fileSize = inFile.tellg();
    if(fileSize == 0)
      throw gError("ApplyVectorFieldFile::setup", "File " + m_fileName + " seems empty!");
    // the modulos (remainders) must be zero
    if(fileSize % sizeof(double) != 0)
      throw gError("ApplyVectorFieldFile::setup", "File " + m_fileName + " does not seem to contain double precision numbers!");
    if(fileSize % (sizeof(double)*3) != 0)
      throw gError("ApplyVectorFieldFile::setup", "File " + m_fileName + " does not seem to contain 3 cartesian components but it must!");
 
    MSG_DEBUG("ApplyVectorFieldFile::setup", "fileSize is " << fileSize);
   
    m_vecField = new double[fileSize/sizeof(double)];
    
    // back to start of file
    inFile.seekg (0, ios::beg);
    inFile.read((char*) m_vecField, fileSize);
    inFile.close();
  }
  else {
    throw gError("ApplyVectorFieldFile::setup", "Can not open file \"" + m_fileName + "\"! Aborting.");
  }

}

void ApplyVectorFieldFile::thermalize(Phase* phase)
{
    FOR_EACH_FREE_PARTICLE_C
      (phase,
       m_colour,
       point_t& vec = __iSLFE->tag.pointByOffset(m_offset);
       point_t& r = __iSLFE->r;
       size_t pSlot = __iSLFE->mySlot;
       vec.x += m_vecField[pSlot*SPACE_DIMS];
       vec.y += m_vecField[pSlot*SPACE_DIMS+1];
       vec.z += m_vecField[pSlot*SPACE_DIMS+2];
    );
}


