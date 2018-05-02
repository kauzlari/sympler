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


#include "pc_free_pcalc.h"

#include "simulation.h"
#include "manager_cell.h"

#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation*) M_PHASE->parent())
#define M_MANAGER M_PHASE->manager()

ParticleCreatorFreePCalc::ParticleCreatorFreePCalc()
  : ParticleCreatorFree()
{
  throw gError("ParticleCreatorFreePCalc::ParticleCreatorFreePCalc()", "Do not use.");
}

ParticleCreatorFreePCalc::ParticleCreatorFreePCalc(Boundary *boundary)
  : ParticleCreatorFree(boundary)
{
  init();
}


ParticleCreatorFreePCalc::~ParticleCreatorFreePCalc()
{
}



//--- Methods ---

void ParticleCreatorFreePCalc::initTransform()
{
  /* Let's hope the size of the box doesn't change anymore after this point. */
  m_box_size = ((Boundary *) m_parent)->boundingBox().size();

  if(m_species == "UNDEF")
  {
    __m_particle.resize(M_MANAGER->nColours());
    for (int i = 0; i < SPACE_DIMS; i++) {
      m_poss[i].resize(M_MANAGER->nColours());
      m_vels[i].resize(M_MANAGER->nColours());
    }
    m_internal_dofs.resize(M_MANAGER->nColours());

    for (size_t c = 0; c < M_MANAGER->nColours(); c++) {
      __m_particle[c].setColour(c);

      for (size_t i = 0; i < SPACE_DIMS; i++) {
        m_vels[i][c] = new ParticleCalculator(&__m_particle[c], c, &m_box_size, m_str_vels[i]);
        m_poss[i][c] = new ParticleCalculator(&__m_particle[c], c, &m_box_size, m_str_poss[i]);
      }

      for (map<string, string>::const_iterator i = m_properties.unknown().begin();
           i != m_properties.unknown().end(); ++i) {
             if (Particle::s_tag_format[c].attrExists(i->first)) {
               dof_info_t info;

               info.offset = Particle::s_tag_format[c].attrByName(i->first).offset;
               info.pc = new ParticleCalculator(&__m_particle[c], c, &m_box_size, i->second);

               m_internal_dofs[c].push_back(info);
             } else {

               string s = string("ParticleCreatorFreePCalc::initTransform") 
		 + string("No internal degree of freedom named '" 
			  + i->first + "' for species '" + M_MANAGER->species(c) + "' exists and no attribute exists. Aborting.");
//                MSG_INFO("ParticleCreatorFreePCalc::initTransform", "No internal degree of freedom named '" + i->first +
//                    "' for species '" + M_MANAGER->species(c) + "' exists and no attribute exists. Aborting.");

               m_properties.throwListIfUnknown(s);	
             }
           }
    }
  }
  else // so, m_species != "UNDEF"
  {
    __m_particle.resize(1);
    for (int i = 0; i < SPACE_DIMS; i++) {
      m_poss[i].resize(1);
      m_vels[i].resize(1);
    }
    m_internal_dofs.resize(1);

    __m_particle[0].setColour(m_colour);

    for (size_t i = 0; i < SPACE_DIMS; i++) {
      m_vels[i][0] = new ParticleCalculator(&__m_particle[0], m_colour, &m_box_size, m_str_vels[i]);
      m_poss[i][0] = new ParticleCalculator(&__m_particle[0], m_colour, &m_box_size, m_str_poss[i]);
    }

    for (map<string, string>::const_iterator i = m_properties.unknown().begin();
         i != m_properties.unknown().end(); ++i) {
           if (Particle::s_tag_format[m_colour].attrExists(i->first)) {
             dof_info_t info;

             info.offset = Particle::s_tag_format[m_colour].attrByName(i->first).offset;
             info.pc = new ParticleCalculator(&__m_particle[0], m_colour, &m_box_size, i->second);
//         info.name = i->first;

             m_internal_dofs[0].push_back(info);
           } else {

             string s = string("ParticleCreatorFreePCalc::initTransform")
	       + string("No internal degree of freedom named '" + i->first +
                 "' for species '" + M_MANAGER->species(m_colour) + "' exists and no attribute exists. Aborting.");
//              MSG_INFO("ParticleCreatorFreePCalc::initTransform", "No internal degree of freedom named '" + i->first +
//                  "' for species '" + M_MANAGER->species(m_colour) + "' exists and no attribute exists. Aborting.");

             m_properties.throwListIfUnknown(s);	
           }
         }
  }
}

void ParticleCreatorFreePCalc::init()
{
  m_properties.setClassName("ParticleCreatorFreePCalc");

  /* Allow unknown properties. Those ones have to be identified later.
  They are used to set the particles degree of freedoms initially. */
  m_properties.allowUnknown();

  for (int i = 0; i < SPACE_DIMS; i++) {
    m_properties.addProperty
        ("pos" + string(1, 'X'+i), PropertyList::STRING, &m_str_poss[i], NULL,
         "Function applied to the " + string(1, 'x'+i) + "-position of each particle "
             "shortly after creation. Example: Can be used to stretch (posX = \"4*posX\") or move "
             "(posX = \"posX+4\") the particles. If the particle is moved out of the given"
             " Boundary, it will NOT be created.");
    m_properties.addProperty
        ("vel" + string(1, 'X'+i), PropertyList::STRING, &m_str_vels[i], NULL,
         "Function applied to the " + string(1, 'x'+i) + "-component of the velocity"
             " of each particle.\nExample1: 'vel" + string(1, 'X'+i) + "=\"5\"' sets the " 
             + string(1, 'x'+i) + "-component of each particle"
             " velocity to 5. This changes the temperature T0 of the particles to"
             " a non-equilibrium temperature T=(2/3)*T0.\nExample2: 'vel" 
             + string(1, 'X'+i) + "=\"vel" + string(1, 'X'+i) 
             + "+5\"' adds a centre of mass velocity of "
             "5 in " + string(1, 'x'+i) + "-direction. The temperature" 
             " is preserved.");

    m_str_poss[i] = "pos" + string(1, 'X'+i);
    m_str_vels[i] = "vel" + string(1, 'X'+i);
  }
}

/* ParticleCalculator */

void ParticleCalculator::AnalyzeId(const string &id)
{
  static map<char, int> space_coords;
	
  if (space_coords.empty()) {
    space_coords['X'] = 0;
    space_coords['Y'] = 1;
    space_coords['Z'] = 2;
  }
	
  if (id == "posX" || id == "posY" || id == "posZ") {
    curr_tok.key = DOUBLE_PTR;
    curr_tok.x_ptr = &m_particle->r[space_coords[id[3]]];
  } else if (id == "velX" || id == "velY" || id == "velZ") {
    curr_tok.key = DOUBLE_PTR;
    curr_tok.x_ptr = &m_particle->v[space_coords[id[3]]];
  } else if (id == "boxX" || id == "boxY" || id == "boxZ") {
    curr_tok.key = DOUBLE_PTR;
    curr_tok.x_ptr = &(*m_box_size)[space_coords[id[3]]];
  } else
    throw gError("ParticleCalculator::AnalyzeId", "Id \"" + id + "\" undefined. " \
        "Possible ids are: posX, posY, posZ, velX, velY, velZ, boxX, boxY, boxZ");
}

void ParticleCreatorFreePCalc::setup()
{
  ParticleCreatorFree::setup();
}
