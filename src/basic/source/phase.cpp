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



#include <unistd.h>
#include <algorithm>

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"
#include "pair_creator.h"
#include "triplet.h"

// #include "valgrind/memcheck.h"

using namespace std;

Phase *DEBUG_g_phase;

Phase* triplet_t::s_phase = NULL;

bool Phase::connectionDeclarationFinished = false;

Phase::Phase(Simulation *simulation)
  : NodeManyChildren((Node*) simulation), 
    // we try to get rid of the following
    // vvSumOld(true), vvSum(0),
     velCMold(true),
     nOfParticles(0), nOfFrozenP(0), m_manager(NULL), m_pairCreator(NULL)
            /*, m_thermostat(NULL)*/
{
   velCM[0] = 0;
   velCM[1] = 0;
   velCM[2] = 0;

  nOfParticles = 0;

  nOfFreeParticlesPerGroup[0] = 0;

//   MSG_DEBUG("Phase::setup", "BEFORE: Valgrindcheck nOfFreeParticlesPerGroup: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(nOfFreeParticlesPerGroup));

  DEBUG_g_phase = this;
  if(triplet_t::s_phase)
    throw gError("Phase::Phase(Simulation *simulation)", "Fatal error: triplet_t::s_phase was already set! Contact the programmers.");
  triplet_t::s_phase = this;

  init();
}


/*virtual*/ Phase::~Phase()
{
//   MSG_DEBUG("Phase::~Phase", "DESTRUCTOR CALLED");
  delete m_manager;
  for(size_t i = 0; i < m_tripletLists.size(); ++i)
    delete m_tripletLists[i];
  //delete m_pairCreator;
}


void Phase::removeParticle(Particle *p)
{
  --nOfParticles;
  --(nOfFreeParticlesPerGroup[p->g]);

  m_particles[p->c].deleteEntry(*p);
}


void Phase::removeParticle(size_t colour, size_t slot)
{
  removeParticle(&m_particles[colour][slot]);
}


Particle *Phase::addParticle(const Particle &particle)
{
  Particle* p;
	p = newParticle(particle.c);
  *p = particle;
	nOfFreeParticlesPerGroup[p->g]++;
  return p;
}


Particle *Phase::addParticle
	(const point_t &pos, const point_t &vel, size_t group, size_t colour)
{
  assert(colour < nColours());

  Particle* p;
	p = newParticle(colour);
  p -> r = pos;
  p -> v = vel;
  p -> g = group;
  p -> setColour(colour);
  nOfFreeParticlesPerGroup[group]++;
  return p;
}


void Phase::addColour(size_t colour)
{
  MSG_DEBUG("Phase::addColour", "called with colour=" << colour 
	    << " and previous nOfColours=" << m_particles.size());
  
  m_particles.resize(colour+1);
  m_frozen_particles.resize(colour+1);
}

// adds particle to end of list and returns interator to added particle
Particle* Phase::newParticle(size_t colour)
{
  Particle &p = m_particles[colour].newEntry();
				
  p.setColour(colour);
  p.clear();
  
  ++nOfParticles;

  return &p;
}

Particle *Phase::addFrozenParticle(const Particle &particle)
{
	Particle &p = m_frozen_particles[particle.c].newEntry();
	p.r = particle.r;
	p.v = particle.v;
	p.g = particle.g;
	p.setColour(particle.c);
	p.isFrozen = 1;

	++nOfFrozenP;

    return &p;
}



Particle *Phase::addFrozenParticle(const point_t &pos, const point_t &vel, size_t group, size_t colour)
{
    Particle &p = m_frozen_particles[colour].newEntry();

    p.r = pos;
    p.v = vel;
    p.g = group;
    p.setColour(colour);
    p.isFrozen = 1;

    ++nOfFrozenP;
//     if(p.mySlot == 1238) MSG_DEBUG("Phase::addFrozenParticle", "NOW 1238");
    return &p;
}

void Phase::setForNewIntegration()
{
  m_manager->clearTags();
}


void Phase::init()
{
  m_properties.setClassName("Phase");

  m_properties.setDescription(
    "The Phase object manages the particle information and initialiazes the simulation."
  );

  BOOLPC(smartCells, m_smartCells, "Specifies whether only those cells should be created, which really lie in the domain. This saves memory but needs much more time during initialisation.");

  BOOLPC(randomPairs, m_randomPairs, "Specifies whether loops over pairs should be performed in a random order.");

  STRINGPC
    (cellFileName, m_cell_filename,
     "Filename of the VTK file holding information about the cell subdivision.");
    
  m_smartCells = false;
  m_randomPairs = false;
  m_cell_filename = "cells.vtk";
  m_particles_assigned = false;
}


void Phase::read(const xmlNode *xmln)
{
  NodeManyChildren::read(xmln);
		
  m_manager = new ManagerCell(this);

  if(!m_boundary)
    throw gError("Phase::read", "No Boundary defined.");
  if(!m_pairCreator)
    throw gError("Phase::read", "No PairCreator defined.");
}


void Phase::setup()
{
  NodeManyChildren::setup();

  if (((Simulation*) m_parent)->randomize())
    m_rng.setSeed(getpid());
  else
    m_rng.setSeed(RNG_DEFAULT_SEED);

//   MSG_DEBUG("Phase::setup", "Valgrindcheck *this: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(*this));

}

void Phase::invalidate()
{
    velCMold = true;
}

// Computes the center of mass velocity
// we compute ONE velCM for all colours
void Phase::computeVelCM()
{
  velCMPerGroup.clear();

  velCM[0] = 0;
  velCM[1] = 0;
  velCM[2] = 0;

  FOR_EACH_FREE_PARTICLE
    (this,
     //        int group = m_particles[*particle].g;

     velCM += __iSLFE->v;

     velCMPerGroup[__iSLFE->g] += __iSLFE->v;
    );

  for (int j = 0; j < SPACE_DIMS; j++) {
    velCM[j] /= nOfParticles;
    
    for (map< size_t, point_t>::iterator g = velCMPerGroup.begin();
	 g != velCMPerGroup.end(); g++)
      g->second[j] /= nOfFreeParticlesPerGroup[g->first];
  }

  velCMold = false;
}


Node *Phase::instantiateChild(const string &name)
{
  Node *node;

  if(Boundary_Factory::exists(name)) {
    node = m_boundary = Boundary_Factory::byName(name).instantiate(this);
  }

  else if (PairCreator_Factory::exists(name)) {
    node = m_pairCreator = PairCreator_Factory::byName(name).instantiate(this);
  } 

  else
    throw gError
      ("Phase::instantiateChild", "'" + name + "' not found.");  
  
  return node;
}


void Phase::assignParticlesToCells()
{
    MSG_DEBUG("Phase::assignParticlesToCells", "Writing cell subdivision information to '"
              << m_cell_filename << "'.");
    m_manager->toVTK(m_cell_filename);

    MSG_DEBUG("Phase::assignParticlesToCells", "Sorting particles to cells.");

    m_manager->clearAll();

    FOR_EACH_FREE_PARTICLE
        (this,
         Cell *cell = m_manager->findCell(__iSLFE->r);
         if (cell) {
           cell->injectFree(c, __iSLFE);
         } else {
             printf("r = (%f, %f, %f)\n"
                    "v = (%f, %f, %f)\n"
                    "g = %lu\nslot = %lu\ncolour = %lu\n",
                    __iSLFE->r.x, __iSLFE->r.y, __iSLFE->r.z,
                    __iSLFE->v.x, __iSLFE->v.y, __iSLFE->v.z,
                    (unsigned long) __iSLFE->g, (unsigned long) __iSLFE->mySlot, (unsigned long) c
                    ); 

            throw gError("Phase::assignParticlesToCells", "No cell for free particle.");
        }
    );

    for (size_t c = 0; c < nColours(); c++)
      FOR_EACH_FROZEN_PARTICLE
        (this, c,
         Cell *cell = m_manager->findCell(__iSLFE->r);
         if (cell) {
           cell->injectFrozen(c, __iSLFE);
         } else {
           printf("r = (%f, %f, %f)\n"
                  "v = (%f, %f, %f)\n"
                  "g = %lu\nslot = %lu\ncolour = %lu\n",
                  __iSLFE->r.x, __iSLFE->r.y, __iSLFE->r.z,
                  __iSLFE->v.x, __iSLFE->v.y, __iSLFE->v.z,
                  (unsigned long) __iSLFE->g, (unsigned long) __iSLFE->mySlot, (unsigned long) c
           ); 
           
           throw gError("Phase::assignParticlesToCells", "No cell for frozen particle.");
         }
    );

    vector<Cell*>::iterator cells_end = m_manager->cells().end();
    for (vector<Cell*>::iterator i = m_manager->cells().begin(); i != cells_end; i++)
        (*i)->commitInjections();    

    countParticles();

    m_particles_assigned = true;

//    m_manager->initActiveLinks();
}



Boundary *Phase::boundary()
{
    return m_boundary;
}

PairCreator *Phase::pairCreator()
{
    return m_pairCreator;
}

void Phase::invalidatePositions(IntegratorPosition *integrator)
{
  m_manager->invalidatePositions(integrator);
  m_pairCreator->invalidatePositions();
}


void Phase::countParticles()
{
  nOfParticles = 0;
  nOfFreeParticlesPerGroup.clear();

  FOR_EACH_FREE_PARTICLE
    (this,
     nOfParticles++;
     //         nOfFreeParticlesPerGroup[m_particles[*particle].g]++;
     nOfFreeParticlesPerGroup[__iSLFE->g]++;
    );

  for(size_t c = 0; c < m_particles.size(); ++c)
  {
    MSG_DEBUG("Phase::countParticles", "number of particles for colour " << c << ": " << m_particles[c].size());
    MSG_DEBUG("Phase::countParticles", "number of frozen particles for colour " << c << ": " << m_frozen_particles[c].size());
  }

//    assert(nOfParticles == m_particles.size());
}


size_t Phase::nColours() const
{
  return m_manager->nColours();
}


void Phase::writeRestartFile(string name)
{
  ofstream pos(name.c_str());
    
  pos.precision(8);
    
  for (size_t c = 0; c < m_manager->nColours(); c++) {
    pos << m_manager->species(c) << " ";
    for (size_t i = 0; i < Particle::s_tag_format[c].rows(); i++) {
      DataFormat::attribute_t attr = Particle::s_tag_format[c].attrByIndex(i);
    
      if (attr.persistent)
	pos << attr.name << " ";
    }
    pos << "!!!" << endl;
  }
  pos << "!!!" << endl;

  FOR_EACH_FREE_PARTICLE
    (this,
     pos << m_manager->species(__iSLFE->c) << " " << "free" << " " << __iSLFE->r.x << " " << __iSLFE->r.y << " " << __iSLFE->r.z 
     << " " << __iSLFE->v.x << " " << __iSLFE->v.y << " " << __iSLFE->v.z;

     for (size_t j = 0; j < Particle::s_tag_format[c].rows(); j++) {
       if (Particle::s_tag_format[c].attrByIndex(j).persistent)
	 pos << " " << __iSLFE->tag.toStringByIndex(j);
     }

     pos << endl;
     );
  
  FOR_EACH_FROZEN_PARTICLE_ALL_C
      (this,
       pos << m_manager->species(__iSLFE->c) << " " << "frozen" << " " << __iSLFE->r.x << " " << __iSLFE->r.y << " " << __iSLFE->r.z 
       << " " << __iSLFE->v.x << " " << __iSLFE->v.y << " " << __iSLFE->v.z;

       for (size_t j = 0; j < Particle::s_tag_format[c].rows(); j++) {
         if (Particle::s_tag_format[c].attrByIndex(j).persistent)
  	 pos << " " << __iSLFE->tag.toStringByIndex(j);
       }

       pos << endl;
       );


  pos << "!!!" << endl;
    
  pos.close();
}


void Phase::addTripletAndWrite(Particle *p1, Particle *p2, Particle *p3, size_t listIndex) {
	triplet_t tmp;
	ofstream file_stream;
	string file_name;
	tmp.a = p1;
	tmp.b = p2;
	tmp.c = p3;
	m_tripletLists[listIndex]->push_back(tmp);
// 	MSG_DEBUG("Phase::addTripletAndWrite", "added triplet (" << p1->mySlot << ", " << p2->mySlot << ", " << p3->mySlot << ")");
	file_name = ((Simulation*) m_parent)->name();
	file_name += "_connector.con";

	file_stream.open(file_name.c_str(), ios_base::app);
	file_stream.seekp(0, std::ios_base::end);

	// first connection in file? Then write the line "!!!" first
	if(!Phase::connectionDeclarationFinished) {
	  Phase::connectionDeclarationFinished = true;
	  file_stream << "!!!" << endl;
	}

	file_stream << "triplet " << tripletListName(listIndex) << " " << p1 -> mySlot << " ";
	if (p1->isFrozen)
		file_stream << "frozen" << " ";
	else
		file_stream << "free" << " ";
	file_stream << p2 -> mySlot << " ";
	if (p2->isFrozen)
		file_stream << "frozen" << " ";
	else
		file_stream << "free" << " ";
	file_stream << p3 -> mySlot << " ";
	if (p3->isFrozen)
		file_stream << "frozen" << " "<< endl;
	else
		file_stream << "free" << " "<< endl;
	file_stream.close();

}

void Phase::addTriplet(Particle *p1, Particle *p2, Particle *p3, size_t listIndex) {
	triplet_t tmp;
	ofstream file_stream;
	string file_name;
	tmp.a = p1;
	tmp.b = p2;
	tmp.c = p3;
	m_tripletLists[listIndex]->push_back(tmp);
// 	MSG_DEBUG("Phase::addTriplet", "added triplet (" << p1->mySlot << ", " << p2->mySlot << ", " << p3->mySlot << ")");
}


/*!
 * Returns the string identifier of the \a tripletList from \a m_tripletLists with index \a listIndex
 */
string Phase::tripletListName(size_t listIndex)
{
  if(listIndex < m_tripletLists.size())
    return m_tripletListNames[listIndex];
  else 
    throw gError("Phase::tripletListName"+FILE_INFO, "Requested index out of bounds! Contact programmer.");
}

/*!
 * Returns the \a tripletList from \a m_tripletLists with string identifier \a name
 */
tripletList* Phase::returnTripletList(string name) {
//   MSG_DEBUG("Phase::returnTripletList", "CALLED");
  size_t i = tripletListIndex(name);
  return m_tripletLists[i];
  }


/*!
 * Returns the index of the \a tripletList from \a m_tripletLists with string identifier \a name
 */
size_t Phase::tripletListIndex(string name)
{
  int tmp = -1;
//   MSG_DEBUG("Phase::tripletListIndex", "Nof tripletLists: " << m_tripletListNames.size() );
  for(size_t i = 0; i < m_tripletListNames.size(); ++i) {
    if(m_tripletListNames[i] == name) {
//       MSG_DEBUG("Phase::tripletListIndex", "now at name " << m_tripletListNames[i] );
      // the same name cannot appear twice, so the following should be OK
      if(tmp == -1) tmp = i;
      else 
	throw gError("Phase::tripletListIndex"+FILE_INFO, "Two triplet lists with same identifier \"" + name + "\" found!");
    }  
  }
  if(tmp == -1) {
    MSG_DEBUG("Phase::tripletListIndex"+FILE_INFO, "before not-found exception");
    throw gError("Phase::tripletListIndex"+FILE_INFO, "No triplet list with identifier \"" + name + "\" found!"); 
  }
  return size_t(tmp);
}

size_t Phase::createTripletListIndexAndWrite(string name)
{
  if(!Phase::connectionDeclarationFinished) {
    int tmp = -1;
    size_t i;
    for(i = 0; i < m_tripletListNames.size(); ++i) {
      if(m_tripletListNames[i] == name) {
	if(tmp == -1) tmp = i;
	else 
	  throw gError("Phase::createTripletListIndexAndWrite"+FILE_INFO, "Two triplet lists with same identifier \"" + name + "\" found! Shouldn't happen. Contact the programmers");
      }  
    }
    if(tmp == -1) {
      // create
      m_tripletListNames.push_back(name);
      m_tripletLists.push_back(new tripletList());

      // write in file
      // FIXME: one of at least four places where we create the same file_name from scratch: redundant and error-prone! Replace by member m_connections_file?
      string file_name = ((Simulation*) m_parent)->name();
      file_name += "_connector.con";
      ofstream file_stream;
      file_stream.open(file_name.c_str(), ios_base::app);
      file_stream.seekp(0, std::ios_base::end);
      file_stream << "triplet " << name << endl;


      return i; // index of the new list
    }
    
    return size_t(tmp); // index of an existing list we found
  }
  else
    throw gError("Phase::createTripletListIndexAndWrite"+FILE_INFO, "Attempt to create a new triplet list eventhough it seems too late for that. Contact programmer.");
}

size_t Phase::createTripletListIndex(string name)
{
  if(!Phase::connectionDeclarationFinished) {
    int tmp = -1;
    size_t i;
    for(i = 0; i < m_tripletListNames.size(); ++i) {
      if(m_tripletListNames[i] == name) {
	if(tmp == -1) tmp = i;
	else 
	  throw gError("Phase::createTripletListIndex"+FILE_INFO, "Two triplet lists with same identifier \"" + name + "\" found! Shouldn't happen. Contact the programmers");
      }  
    }
    if(tmp == -1) {
      // create
      m_tripletListNames.push_back(name);
      m_tripletLists.push_back(new tripletList());

      return i; // index of the new list
    }
    
    return size_t(tmp); // index of an existing list we found
  }
  else
    throw gError("Phase::createTripletListIndex"+FILE_INFO, "Attempt to create a new triplet list eventhough it seems too late for that. Contact programmer.");
}

/*!
 * Returns the index of the \a tripletList from \a m_tripletLists with string identifier \a name
 */
int Phase::searchTripletListIndex(string name)
{
  int tmp = -1;
  for(size_t i = 0; i < m_tripletListNames.size(); ++i) {
    if(m_tripletListNames[i] == name) {
      if(tmp == -1) tmp = i;
      else throw gError("Phase::tripletListIndex"+FILE_INFO, "Two triplet lists with same identifier \"" + name + "\" found!");
    }
  }
  return tmp;
}

// this function assumes that the Valcalculator does all the registering work itself
void Phase::registerBondedCalc(TripletCalculator* vC)
{
  // this should be a safe check whether this function should be used
  if(!vC->parent())
    throw gError("Phase::registerBondedCalc(ValCalculator* vC)", "The function shouldn't have been called because vC does not belong to a node hierarchy. Contact the programmer.");
  
  m_bondedTripletCalculators_flat.push_back(vC);

//   MSG_DEBUG("Phase::registerBondedCalc", "m_bondedTripletCalculators_flat.size now = " << m_bondedTripletCalculators_flat.size() << ", name[last] = " << m_bondedTripletCalculators_flat[m_bondedTripletCalculators_flat.size()-1]->className());  
}

// this function assumes that the Valcalculator does all the registering work itself
void Phase::registerBondedCalc_0(TripletCalculator* vC)
{
  // this should be a safe check whether this function should be used
  if(!vC->parent())
    throw gError("Phase::registerBondedCalc_0(ValCalculator* vC)", "The function shouldn't have been called because vC does not belong to a node hierarchy. Contact the programmer.");
  
  m_bondedTripletCalculators_flat_0.push_back(vC);

//   MSG_DEBUG("Phase::registerBondedCalc_0", "m_bondedTripletCalculators_flat_0.size now = " << m_bondedTripletCalculators_flat_0.size() << ", name[last] = " << m_bondedTripletCalculators_flat_0[m_bondedTripletCalculators_flat_0.size()-1]->className());  
}


bool Phase::findStages()
{    
  bool finished = true;

  for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculators_flat.begin(); i != m_bondedTripletCalculators_flat.end(); ++i)
    // must be findStage() first so that it is definitely executed
    finished = (*i)->findStage() && finished; 

  return finished;
}

bool Phase::findStages_0()
{    
  bool finished = true;

  for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculators_flat_0.begin(); i != m_bondedTripletCalculators_flat_0.end(); ++i)
    // must be this order so that findSTage() is definitely executed
    finished = (*i)->findStage_0() && finished; 

  return finished;
}

void Phase::sortStages()
{
  // bonded Calculators
  m_maxBondedStage.resize(m_tripletLists.size());
  m_bondedTripletCalculators.resize(m_tripletLists.size());
  for(size_t icl = 0; icl < m_maxBondedStage.size(); ++icl) 
    m_maxBondedStage[icl] = 0;

  for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculators_flat.begin(); i != m_bondedTripletCalculators_flat.end(); ++i)
    {
      // preliminary cast as long as we don't have the complete hierarchy
      string listName = ((TripletCalculator*) (*i))->listName();
      size_t listIndex = tripletListIndex(listName);

//       MSG_DEBUG("Phase::sortStages", "now bonded Calculator " << (*i)->className() << " with stage = " << (*i)->stage()) << ", for list \"" << listName << "\".";

      int stage = (*i)->stage();
      assert(stage > -1);
      if((size_t) stage > m_maxBondedStage[listIndex]) m_maxBondedStage[listIndex] = stage;
      
      if(m_bondedTripletCalculators[listIndex].size() <= (size_t) stage)
	{
// 	  MSG_DEBUG("Phase::sortStages", "increasing the number of bonded stage-slots to " << stage+1 << " because of " << (*i)->className());

	  m_bondedTripletCalculators[listIndex].resize(stage+1);
	}
#if 0 // this is when we will have a separation of Calculators
#ifdef _OPENMP
      if((m_valCalculatorParts.size() <= (size_t) stage))
	{
	  MSG_DEBUG("Phase::sortStages", "increasing the number of stage-slots to " << stage+1 << " because of " << (*i)->className());
	  m_valCalculatorParts.resize(stage+1);
	}
#endif
#endif
      
#ifndef _OPENMP
      assert((size_t)stage < m_bondedTripletCalculators[listIndex].size());
      m_bondedTripletCalculators[listIndex][stage].push_back(*i);
#else
      // The commented out code will be needed when we introduce a separation of Calculators
      // 	if ((*i)->particleCalculator()) {
      // 	  assert((size_t)stage < m_valCalculatorParts.size());
// 	  m_valCalculatorParts[stage].push_back((ValCalculatorPart*)(*i));  
// 	}
// 	else {
      assert((size_t)stage < m_bondedTripletCalculators[listIndex].size());
      m_bondedTripletCalculators[listIndex][stage].push_back(*i);    
      // 	}
#endif
    } // end of loop over m_bondedTripletCalculators_flat

  for(size_t listIndex = 0; listIndex < m_bondedTripletCalculators.size(); ++listIndex) {
    if(!m_bondedTripletCalculators[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculators[listIndex][s].begin(); i != m_bondedTripletCalculators[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("Phase::sortStages", "found stage " << s << " for " << (*i)->className());
	      }
	  }
      }
    else m_bondedTripletCalculators[listIndex].resize(1);
#if 0 // next will be needed when separating the bonded Calculators
#ifdef _OPENMP
    if(!m_bondedalCalculatorParts[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculatorParts[listIndex][s].begin(); i != m_bondedTripletCalculatorParts[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("Phase::sortStages", "found stage " << s << " for " << (*i)->className());
	      }
	  }
      }
    else m_bondedTripletCalculatorParts[listIndex].resize(1);
#endif
#endif
  }
    
}


void Phase::sortStages_0()
{

  // bonded Calculators
  m_maxBondedStage_0.resize(m_tripletLists.size());
  m_bondedTripletCalculators_0.resize(m_tripletLists.size());
  for(size_t icl = 0; icl < m_maxBondedStage_0.size(); ++icl) 
    m_maxBondedStage_0[icl] = 0;

  for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculators_flat_0.begin(); i != m_bondedTripletCalculators_flat_0.end(); ++i)
    {
      // preliminary cast as long as we don't have the complete hierarchy
      string listName = ((TripletCalculator*) (*i))->listName();
      size_t listIndex = tripletListIndex(listName);

      int stage = (*i)->stage();
      assert(stage > -1);
      if((size_t) stage > m_maxBondedStage_0[listIndex]) m_maxBondedStage_0[listIndex] = stage;
      
      if(m_bondedTripletCalculators_0[listIndex].size() <= (size_t) stage)
	{
	  m_bondedTripletCalculators_0[listIndex].resize(stage+1);
	}
#if 0 // this is when we will have a separation of Calculators
#ifdef _OPENMP
      if((m_valCalculatorParts_0.size() <= (size_t) stage))
	{
	  m_valCalculatorParts_0.resize(stage+1);
	}
#endif
#endif
      
#ifndef _OPENMP
      assert((size_t)stage < m_bondedTripletCalculators_0[listIndex].size());
      m_bondedTripletCalculators_0[listIndex][stage].push_back(*i);
#else
      // The commented out code will be needed when we introduce a separation of Calculators
      // 	if ((*i)->particleCalculator()) {
      // 	  assert((size_t)stage < m_valCalculatorParts_0.size());
      // 	  m_valCalculatorParts_0[stage].push_back((ValCalculatorPart*)(*i));  
      // 	}
      // 	else {
      assert((size_t)stage < m_bondedTripletCalculators_0[listIndex].size());
      m_bondedTripletCalculators_0[listIndex][stage].push_back(*i);    
      // 	}
#endif
    } // end of loop over m_bondedTripletCalculators_flat

  for(size_t listIndex = 0; listIndex < m_bondedTripletCalculators_0.size(); ++listIndex) {
    if(!m_bondedTripletCalculators_0[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage_0[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculators_0[listIndex][s].begin(); i != m_bondedTripletCalculators_0[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("Phase::sortStages_0", "found stage " << s << " for " << (*i)->className());
	      }
	  }
      }
    else m_bondedTripletCalculators_0[listIndex].resize(1);
  
#if 0 // next will be needed when separating the bonded Calculators
#ifdef _OPENMP
    if(!m_bondedalCalculatorParts_0[listIndex].empty())
      {
	for(size_t s = 0; s <= m_maxBondedStage_0[listIndex] /*s_maxStage*/; ++s)
	  {
	    for(vector<TripletCalculator*>::iterator i = m_bondedTripletCalculatorParts_0[listIndex][s].begin(); i != m_bondedTripletCalculatorParts_0[listIndex][s].end(); ++i)
	      {
		MSG_DEBUG("Phase::sortStages_0", "found stage " << s << " for " << (*i)->className());
	      }
	  }
      }
    else m_bondedTripletCalculatorParts_0[listIndex].resize(1);
#endif
#endif
  }

}
