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



#include "meter_pos_vel.h"

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"
#include "pair_creator.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterPosVel> meter_pos_vel("MeterPosVel");

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

#define IDX_TIME 0
#define IDX_DATA_FORMAT 1
#define IDX_POSITIONS 2
#define IDX_VELOCITIES 3
#define IDX_FORCE1 4
#define IDX_FORCE2 5
#define IDX_COLOURS 6


//---- Constructors/Destructor ----

MeterPosVel::MeterPosVel(Simulation *simulation)
  : MeterGrouped(simulation), m_colour(0)
{
    init();
}


/*virtual*/ MeterPosVel::~MeterPosVel()
{
}



//---- Methods ----

/*virtual*/ void MeterPosVel::measureNow(const double& time)
{
  Phase *phase = ((Simulation *) m_parent)->phase();
  //     int dim = phase->returnDim();
  //     double subValue[SPACE_DIMS];
  
  vector_point_sp velocities;
  vector_point_sp positions;
  vector_point_sp force1;
  vector_point_sp force2;
  vector_int_sp colours;
  vector<vector_double_sp> add_doubles;
  vector<vector_point_sp> add_points;
  vector<vector_tensor_sp> add_tensors;

// MSG_DEBUG("MeterPosVel::measureNow", "Writing position file at time " << time);

//   phase->pairCreator()->createDistances();
   

  velocities.alloc();
  positions.alloc();
  colours.alloc();
  force1.alloc();
  force2.alloc();

  if (m_colour == ALL_COLOURS)
    add_doubles.resize(m_add_doubles[0].size());
  else {
    add_doubles.resize(m_add_doubles[m_colour].size());
    add_points.resize(m_add_points[m_colour].size());
    add_tensors.resize(m_add_tensors[m_colour].size());
  }

  for (vector<vector_double_sp>::iterator i = add_doubles.begin(); i != add_doubles.end(); i++) {
    i->alloc();
  }
  for (vector<vector_point_sp>::iterator i = add_points.begin(); i != add_points.end(); i++) {
    i->alloc();
  }
  for (vector<vector_tensor_sp>::iterator i = add_tensors.begin(); i != add_tensors.end(); i++) {
    i->alloc();
  }

  point_t null_vector = { { { 0, 0, 0 } } };
  tensor_t null_tensor; /* Is initialized to zero. */

  if(m_withFrozen)
  {
    FOR_EACH_PARTICLE_C
      (
        phase, m_colour,

        positions->push_back(__iSLFE->r);
        velocities->push_back(__iSLFE->v);
        force1->push_back(__iSLFE->force[0]);
        force2->push_back(__iSLFE->force[1]);
        colours->push_back(__iSLFE->c);
     
        for (size_t j = 0; j < m_add_doubles[c].size(); ++j)
          if (m_add_doubles[c][j] < 0)
            add_doubles[j]->push_back(0);
          else
            add_doubles[j]->push_back(__iSLFE->tag.doubleByIndex(m_add_doubles[c][j]));

        for (size_t j = 0; j < m_add_points[c].size(); ++j)
          if (m_add_points[c][j] < 0)
            add_points[j]->push_back(null_vector);
          else
            add_points[j]->push_back(__iSLFE->tag.pointByIndex(m_add_points[c][j]));

        for (size_t j = 0; j < m_add_tensors[c].size(); ++j)
          if (m_add_tensors[c][j] < 0)
            add_tensors[j]->push_back(null_tensor);
          else
            add_tensors[j]->push_back(__iSLFE->tag.tensorByIndex(m_add_tensors[c][j]));
      );
  
  }
  else
  {
    FOR_EACH_FREE_PARTICLE_C
      (
        phase, m_colour,
         
        positions->push_back(__iSLFE->r);
        velocities->push_back(__iSLFE->v);
        force1->push_back(__iSLFE->force[0]);
        force2->push_back(__iSLFE->force[1]);
        colours->push_back(__iSLFE->c);

        for (size_t j = 0; j < m_add_doubles[c].size(); ++j)
          if (m_add_doubles[c][j] < 0)
            add_doubles[j]->push_back(0);
          else
            add_doubles[j]->push_back(__iSLFE->tag.doubleByIndex(m_add_doubles[c][j]));

        for (size_t j = 0; j < m_add_points[c].size(); ++j)
          if (m_add_points[c][j] < 0)
            add_points[j]->push_back(null_vector);
          else
          {
            add_points[j]->push_back(__iSLFE->tag.pointByIndex(m_add_points[c][j]));
          }
        for (size_t j = 0; j < m_add_tensors[c].size(); ++j)
        {
          if (m_add_tensors[c][j] < 0)
          {
            add_tensors[j]->push_back(null_tensor);
          }
          else
          {
            add_tensors[j]->push_back(__iSLFE->tag.tensorByIndex(m_add_tensors[c][j]));
          }
        }
      );
  }
    
  data_sp data = m_format.newData();
//   MSG_DEBUG("MeterPosVel::measureNow", "m_format.rows()=" << m_format.rows());
  data->doubleByIndex(IDX_TIME) = time;
  data->intByIndex(IDX_DATA_FORMAT) = DF_POINTS;
  data->vectorPointByIndex(IDX_POSITIONS) = positions;
  data->vectorPointByIndex(IDX_VELOCITIES) = velocities;
  data->vectorPointByIndex(IDX_FORCE1) = force1;
  data->vectorPointByIndex(IDX_FORCE2) = force2;
  data->vectorIntByIndex(IDX_COLOURS) = colours;

  size_t i;
  for (i = 0; i < add_doubles.size(); ++i) {
    data->vectorDoubleByIndex(IDX_COLOURS+1+i) = add_doubles[i];
  }

  size_t j;
  for (j = 0; j < add_points.size(); ++j) {
    data->vectorPointByIndex(IDX_COLOURS+1+i+j) = add_points[j];
  }

  i += j;
  for (j = 0; j < add_tensors.size(); ++j) {
    data->vectorTensorByIndex(IDX_COLOURS+1+i+j) = add_tensors[j];
/*    MSG_DEBUG("MeterPosVel::measureNow", "adding tensor to data: i=" << i << ", j=" << j  << ", first tensor = " << (*(add_tensors[j]))[0]);*/
    
  }
  
  distribute(data);
}


void MeterPosVel::init()
{
  m_properties.setClassName("MeterPosVel");

  m_properties.setDescription(
    "Outputs the positions and velocities of all particles "
    "of species 'species'. For additional complete output of all user-defined particle-data, each species needs a separate MeterPosVel."
    "Output "
    "columns are then 'time', 'positions', 'velocities', 'colour' and the additional" 
    " attributes of the species. If the attribute 'species' remains undefined, all species are outputted, but WITHOUT user-defined particle-data."
  );

  STRINGPC
    (species,
     m_species,
     "Defines the species, which should be tracked.");
    
  BOOLPC
    (withFrozen,
     m_withFrozen,
     "Specifies, whether frozen particles (usually wall particles), should be saved too.");
    
  m_species = "UNDEF";    
  m_withFrozen = false;
    
/*  MSG_DEBUG("MeterPosVel::init", "before addAttribute: m_format.rows()=" 
    << m_format.rows());*/
  
  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute(IDDF_DATA_FORMAT, DataFormat::INT);
  m_format.addAttribute("position", DataFormat::VECTOR_POINT);
  m_format.addAttribute("velocity", DataFormat::VECTOR_POINT);
  m_format.addAttribute("force1", DataFormat::VECTOR_POINT);
  m_format.addAttribute("force2", DataFormat::VECTOR_POINT);
  m_format.addAttribute("colour", DataFormat::VECTOR_INT);
}


void MeterPosVel::setup()
{
  Meter::setup();
}

void MeterPosVel::aboutToStart()
{
  static char *possible_attributes[] = { "local_density", "internal_energy", NULL };

  m_add_doubles.resize(M_MANAGER->nColours());
  m_add_points.resize(M_MANAGER->nColours());
  m_add_tensors.resize(M_MANAGER->nColours());

  if (m_species == "UNDEF")
    m_colour = ALL_COLOURS;
  else
    m_colour = M_MANAGER->getColour(m_species);

  /* If we output all colours, look for the predefined attributes given in
     possible_attributes. */
  if (m_colour == ALL_COLOURS) {
    for (char **cura = possible_attributes; *cura; cura++) {
      m_format.addAttribute(*cura, DataFormat::VECTOR_DOUBLE);

      for (size_t c = 0; c < M_MANAGER->nColours(); ++c) {
        if (Particle::s_tag_format[c].attrExists(*cura)) {
          m_add_doubles[c].push_back(Particle::s_tag_format[c].attrByName(*cura).index);
        } else {
          m_add_doubles[c].push_back(-1);
        }
      }
    }
  } else {
    m_add_doubles[m_colour].clear();
    m_add_points[m_colour].clear();
    m_add_tensors[m_colour].clear();

    for (size_t i = 0; i < Particle::s_tag_format[m_colour].rows(); ++i) {
      DataFormat::attribute_t attr;

      attr = Particle::s_tag_format[m_colour].attrByIndex(i);

      /* For now, we only support outputting doubles, points and tensors. */
      if (attr.datatype == DataFormat::DOUBLE) {
        m_add_doubles[m_colour].push_back(i);
      } else if (attr.datatype == DataFormat::POINT) {
        m_add_points[m_colour].push_back(i);
      } else if (attr.datatype == DataFormat::TENSOR) {
        m_add_tensors[m_colour].push_back(i);
      }
    }

    for (size_t i = 0; i < m_add_doubles[m_colour].size(); ++i) {
      DataFormat::attribute_t attr;

      attr = Particle::s_tag_format[m_colour].attrByIndex(m_add_doubles[m_colour][i]);

      m_format.addAttribute(attr.name, DataFormat::VECTOR_DOUBLE);
    }

    for (size_t i = 0; i < m_add_points[m_colour].size(); ++i) {
      DataFormat::attribute_t attr;

      attr = Particle::s_tag_format[m_colour].attrByIndex(m_add_points[m_colour][i]);

      m_format.addAttribute(attr.name, DataFormat::VECTOR_POINT);
    }

    for (size_t i = 0; i < m_add_tensors[m_colour].size(); ++i) {
      DataFormat::attribute_t attr;

      attr = Particle::s_tag_format[m_colour].attrByIndex(m_add_tensors[m_colour][i]);

      m_format.addAttribute(attr.name, DataFormat::VECTOR_TENSOR);
    }
  }

}

