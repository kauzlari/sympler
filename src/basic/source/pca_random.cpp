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



#include "particle_cache.h"
#include "pca_random.h"
#include "simulation.h"
#include "manager_cell.h"



const SymbolRegister<PCaRandom> pca_random("PCaRandom");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


PCaRandom::PCaRandom(/*Node*/Simulation* parent/*size_t colour*/)
/*: m_colour(colour), m_stage(0)*/
  : ParticleCache(parent)
{
  m_stage=1; // no stage dependency just set for no reason to 1
  m_rng_type = gsl_rng_default;
  m_rgen = gsl_rng_alloc (m_rng_type);
  m_datatype = DataFormat::POINT;
  init();
}

PCaRandom::~PCaRandom()
{
}

void PCaRandom::setup() {
  MSG_DEBUG("PCaRandom::setup","In setup");
  
  m_colour = M_MANAGER->getColour(m_species);

  MSG_DEBUG("PCaRandom::setup","Using distribution: " + m_distribution);
  MSG_DEBUG("PCaRandom::setup","Using data type: + " + m_type);
  MSG_DEBUG("PCaRandom::setup","Using species: + " + m_species);
  MSG_DEBUG("PCaRandom::setup","Using symbol name: + " + m_symbolName);
  MSG_DEBUG("PCaRandom::setup","Using colour number name: + " + m_colour);
  
  if (m_species=="default") 
    throw gError("PCaRandom::setup","You must assign a species.");
  if (!(m_distribution=="gauss"))
    throw gError("PCaRandom::setup","Unsupported distribution:" + m_distribution);
  //if ( !((m_type=="scalar") || (m_type=="vector") || (m_type=="tensor")))
  if (!(m_type=="tensor"))         
    throw gError("PCaRandom::setup","Data type is not supported:" + m_type + " Currently supported: tensor");
  if (m_symbolName=="default")
    throw gError("PCaRandom::setup","You must assign a valid and unique name for 'symbol'.");
  if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
    throw gError("PCaRandom::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed in Eigensystem.");

  Particle::registerCache(this);
}

void PCaRandom::init()
{
  MSG_DEBUG("PCaRandom::init","In init");

  m_properties.setClassName("PCaRandom");
  m_properties.setName("PCaRandom");
  m_properties.setDescription("Returns gaussian etc distributed random values");

  m_seed=0;
  m_sigma=1.0;
  m_distribution="undefined";
  m_type="undefined";
  m_symbolName="undefined";
  
  STRINGPC(distribution, m_distribution, "Specifies which random number distribution to use.");
  DOUBLEPC(sigma, m_sigma, 0, "Set standard deviation for the gauss random number generator.");
  STRINGPC(type, m_type, "Set the type of the random number: scalar, vector or tensor.");
  INTPC(seed, m_seed, 0, "Specify a seed for the random generator. If 0, the system time will be used.");
  
  STRINGPC
    (symbol, m_symbolName,
     "Symbol name for the random object.");
  

 
  MSG_DEBUG("PCaRandom::init","in setup");
  MSG_DEBUG("PCaRandom::init","Using distribution: " + m_distribution);
  MSG_DEBUG("PCaRandom::init","Using data type: + " + m_type);
  MSG_DEBUG("PCaRandom::init","Using spceis: + " + m_species);
  MSG_DEBUG("PCaRandom::init","Using symbol name: + " + m_symbolName);
}


void PCaRandom::computeCacheFor(Particle* p) {
  //todo for vector retval is simple ref
  //for other: smartpointer ?
  if (m_type=="scalar") {
    //    MSG_DEBUG("PCaRandom::computeCache","Calc. scalar property.");
    double& value = p->tag.doubleByOffset(m_offset);
    value = gsl_ran_gaussian (m_rgen, m_sigma);
  } 
  else if (m_type=="vector") {
    //MSG_DEBUG("PCaRandom::computeCache","Calc. vector property.");
    vector_double_sp& value = p->tag.vectorDoubleByOffset(m_offset);
    for(size_t i = 0; i < SPACE_DIMS; ++i) {
      //todo
    }
  } 
  else if (m_type=="tensor") {
    //MSG_DEBUG("PCaRandom::computeCache","Calc. tensor property.");
    tensor_t& tensor = p->tag.tensorByOffset(m_offset);
    for(size_t i = 0; i < SPACE_DIMS; ++i)
      for(size_t j = 0; j < SPACE_DIMS; ++j)
	tensor(i,j)=gsl_ran_gaussian (m_rgen, m_sigma);
  }
}


void PCaRandom::registerWithParticle() {
  MSG_DEBUG("PCaRandom","Registering additional DOF.");
  if (m_type=="scalar") {
    MSG_DEBUG("PCaRandom::registerWithParticle","Reg. scalar property.");
    m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, DataFormat::DOUBLE, false, m_symbolName).offset;
  } else if (m_type=="vector") {
    MSG_DEBUG("PCaRandom::registerWithParticle","Reg. vector property.");
    m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, DataFormat::VECTOR_DOUBLE, false, m_symbolName).offset;
  } else if (m_type=="tensor") {
    MSG_DEBUG("PCaRandom::registerWithParticle","Reg. tensor property.");
    m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, DataFormat::TENSOR, false, m_symbolName).offset;
  } else {
    throw new gError("Unsupported type:" + m_type);
  }
}

