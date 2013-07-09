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



#include "controller.h"

#include "function.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include <ctime>
#include "particle_cache.h"
#include "threads.h"
#include "triplet.h"

#ifdef _OPENMP
  #include "omp.h"
#endif


using namespace std;
/*double*/int Controller::time_for_parallel;
/*double*/int Controller::time_for_parallel1;
/*double*/int Controller::time_for_parallel2;
/* Report every 10 seconds. */
#define REPORT_TIME 10

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()


Controller::Controller(Simulation *simulation): NodeManyChildren((Node*) simulation), m_t(0), m_force_index(0)
{
  init();
}

Controller::~Controller()
{
}


/*--- Methods ---*/

void Controller::init()
{
  m_properties.setClassName("Controller");

  m_properties.setDescription(
    "Controller controls the time integration. Its children "
    "are Integrators."
  );

  INTPC
      (timesteps, m_timesteps, 0,
       "Number of timesteps for the simulation to run.");
  INTPC
      (statusEvery, m_statusEvery, 0,
       "How often (in timesteps) should the progress of the simulation be printed?");
  DOUBLEPC
    (dt, m_dt, 0,
     "Time increment. Total simulation time is dt*timesteps.");

  m_timesteps = 100;
  m_statusEvery = 100;
  m_dt = 0.005;

  m_pairCreateTime = 0;
  m_pairForceTime = 0;
  m_otherForceTime = 0;
  m_integrateTime = 0;
  m_initTime = 0;
  m_partSymbolsTime = 0;
  m_bondedSymbolsTime = 0;
  m_nonBondedSymbolsTime = 0;
  m_SpecialTime = 0;
}


void Controller::run() { 
  //      MSG_DEBUG("Controller::run", M_SIMULATION->phase()->manager()->firstLink()->next->actsOn().first << "actsOnInfo starting");

  clock_t tInitStart, tInitEnd;
//   time_t tInitStart, tInitEnd;
  tInitStart = clock();
//   std::time(&tInitStart);
  
  Phase* phase = ((Simulation*) m_parent)->phase();
  
  MSG_DEBUG("Controller:run", "Hard compiling all functions.");
  
  // compile all the functions from the input file
  for(set<Function*>::iterator funcs = Function::toCompile.begin();
      funcs != Function::toCompile.end(); ++funcs)
    (*funcs) -> compile();
  // don't need that anymore
  Function::toCompile.clear();
  MSG_DEBUG("Controller:run", "Hard compiling DONE.");

  double timestep;//, minsize;

  // FIXME: 7 calls via M_SIMULATION follow: make one M_SIMULATION-call out of that
  M_SIMULATION->setSymbolStages();

  // this may be called before setupAfterParticleCreation() because the (empty) connected lists have already been created which is important for setting the stages of the tripletLists
  M_SIMULATION->sortSymbolStages();

#ifdef _OPENMP
  // next to be called AFTER sortSymbolStages()
  // next to be called BEFORE createParticles()
  M_SIMULATION->setupCopyVectors();
#endif

  M_SIMULATION->phase()->boundary()->setup(M_SIMULATION, M_SIMULATION->phase()->manager());
  
  // next to be called AFTER setupCopyVectors()
  // next to be called BEFORE setupAfterParticleCreation()
  M_SIMULATION->phase()->boundary()->createParticles();
  MSG_DEBUG("Controller::run", "all particles created");

  M_SIMULATION->phase()->assignParticlesToCells();
  //   MSG_DEBUG("Controller::run", "assignParticlesToCells() finished");

  // moved here (old place see below, 2011-03-17)
  // FIXME: probably we can merge this with setupAfterParticleCreation(), but the reason I moved this before setupAfterParticleCreation() was that I need the setup of the Meters BEFORE the setup of its Postprocessors. So, if this is still fulfilled, THEN AND ONLY THEN you can merge!!!
  M_SIMULATION->setupMeters();

  // next to be called AFTER createParticles()
  FOR_EACH
    (list<Node*>,
     m_toSetupAfterParticleCreation,
     (*__iFE)->setupAfterParticleCreation();
     );
  
  m_toSetupAfterParticleCreation.clear();

  // a consistency check for Symbols; currently (2009-08-10) 
  // just needed for tripletCalculators
  // FIXME: again, the need shows up to generalise things like 
  // setupAfterParticleCreation() with a sort of dependency-tree 
  // and recursive setup calls
  for(vector<Symbol*>::iterator symIt =  M_SIMULATION->symbols()->begin(); symIt != M_SIMULATION->symbols()->end(); ++symIt) {
    (*symIt)->checkConsistency();
  }

  MSG_DEBUG("Controller::run", "Showing information stored for each particle:");

  for(size_t c = 0; c < M_MANAGER->nColours(); ++c)
    {
      cout << "#####" << " species " << M_SIMULATION->phase()->manager()->species(c) << " ##############################################" << endl;
      for (size_t i = 0; i < Particle::s_tag_format[c].rows(); i++) {
	cout << Particle::s_tag_format[c].attrByIndex(i).name << endl;
      }
      cout << "##############################################################" << endl;
    }
  
  MSG_DEBUG("Controller::run", "Showing information stored for each pair:");
  
  FOR_EACH_COLOUR_PAIR
    (
     M_SIMULATION->phase()->manager(),
     cout << "#####" << " CP(" << cp->firstSpecies() << ", " << cp->secondSpecies()
     << ") ###############################################" << endl;
     for (size_t i = 0; i < cp->tagFormat().rows(); i++) {
       cout << cp->tagFormat().attrByIndex(i).name << endl;
     }
     cout << "##############################################################" << endl;
     );

  FOR_EACH_FREE_PARTICLE
    (M_SIMULATION->phase(),
     for (int j = 0; j < FORCE_HIST_SIZE; j++)
       __iSLFE->force[j].assign(0);
     __iSLFE->clear(0);
     );

  FOR_EACH
    (list<Node*>,
     m_children,
     ((Integrator*) *__iFE)->isAboutToStart();
     );
  /*    computation of derived quantities: currently (2005/04/22) only necessary for 
	computation of temperature in IntegratorEnergy;
	the temperature (and other quantities) may be necessary for computation of forces,
	that's why we call it here
  */
  FOR_EACH
    (
     list<Node*>, m_children,
     ((Integrator*) *__iFE)->deriveQuantities();
     );
  
  M_SIMULATION->phase()->setForNewIntegration();

  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_initTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_initTime += difftime(tInitEnd, tInitStart);

  tInitStart = clock();
//   std::time(&tInitStart);

  FOR_EACH_COLOUR_PAIR
    ( 
     phase->manager(),
     cp->updateConnectedDistances();    
     );

  M_PAIRCREATOR->createDistances();

  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_pairCreateTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_initTime += difftime(tInitEnd, tInitStart);

  runSymbols();

  tInitStart = clock();
//   std::time(&tInitStart);


#pragma omp parallel for
  for (int t = 0; t < global::n_threads; ++t) {
    vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
    for(vector<ColourPair*>::iterator __cp = M_MANAGER->colourPairs().begin(); __cp != __end; ++__cp) {                                                
      ColourPair *cp = *__cp; 

      vector<GenF*>::iterator pairForcesBegin = cp->pairForces()->begin();
      vector<GenF*>::iterator pairForcesEnd = cp->pairForces()->end();
      
      if(cp->freePairsRandom(t).size())	{
	    for (PrimitiveSLEntry<size_t> *i = cp->freePairsRandom(t).first(); i != NULL; i = i->next) {
	      Pairdist* pair = &(cp->freePairs()[t][i->m_val]);
	    
	      for (vector<GenF*>::iterator force = pairForcesBegin; force != pairForcesEnd; ++force) {
	        (*force)->invalidate();
#ifndef _OPENMP
	        (*force)->computeForces(pair, m_force_index);
#else
	        (*force)->computeForces(pair, m_force_index, t);
#endif
	      }
	    }
	  }
      else {
	    for (Pairdist *i = cp->freePairs()[t].first(); i != NULL; i = i->next) {
	    
	      Pairdist* pair = i;
	    
	      for (vector<GenF*>::iterator force = pairForcesBegin; force != pairForcesEnd; ++force) {
	        (*force)->invalidate();
#ifndef _OPENMP
	        (*force)->computeForces(pair, m_force_index);
#else
	        (*force)->computeForces(pair, m_force_index, t);
#endif
	      
	      }
	    } 
      }
      
      if(cp->frozenPairsRandom(t).size()) {
	    for (PrimitiveSLEntry<size_t> *i = cp->frozenPairsRandom(t).first(); i != NULL; i = i->next) {
	      Pairdist* pair = &(cp->frozenPairs()[t][i->m_val]);
	    
	      for (vector<GenF*>::iterator force = pairForcesBegin; force != pairForcesEnd; ++force) {
	        (*force)->invalidate();
#ifndef _OPENMP
	        (*force)->computeForces(pair, m_force_index);
#else
	        (*force)->computeForces(pair, m_force_index, t);
#endif
	      }
	    
	    }
	  }
      else {
	    for (Pairdist *i = cp->frozenPairs()[t].first(); i != NULL; i = i->next) {
	      Pairdist* pair = i;
	    
	      for (vector<GenF*>::iterator force = pairForcesBegin; force != pairForcesEnd; ++force) {
	        (*force)->invalidate();
#ifndef _OPENMP
	        (*force)->computeForces(pair, m_force_index);
#else
	        (*force)->computeForces(pair, m_force_index, t);
#endif
	      }
	    }
	  }  
    }
  } // end of for (int t = 0; t < global::n_threads; ++t)

#ifdef _OPENMP
  for (int t = 0; t < global::n_threads; t++) {
    list<Node*>::iterator i_begin = m_children.begin();
    list<Node*>::iterator i_end = m_children.end();
    for(list<Node*>::iterator integr = i_begin; integr != i_end; ++integr) {
      size_t _c = ((Integrator*)(*integr))->colour();
      
      FOR_EACH_PARTICLE_C 
        (M_SIMULATION->phase(), _c,
	 ((Integrator*)(*integr))->mergeCopies(__iSLFE, t, m_force_index);
	 );      
    }
  }

#endif
  
  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_pairForceTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;

  tInitStart = clock();
//   std::time(&tInitStart);

  
  for(size_t c = 0; c < M_PHASE->nColours(); ++c) {
    for (Particle *i = M_PHASE->particles(c).first(); i != NULL; i = i->next) {
      for (vector<GenF*>::iterator force = (*M_SIMULATION->particleForces())[c]->begin(); force != (*M_SIMULATION->particleForces())[c]->end(); ++force) {
	    (*force)->invalidate();
	    (*force)->computeForces(i, m_force_index);
	//	MSG_DEBUG("Controller::run()", "(*force)->name() = " << (*M_SIMULATION->particleForces())[c]->size());
      }
    }
  }

  FOR_EACH
    (vector<GenF*>,
     (*M_SIMULATION->otherForces()),
     (*__iFE)->invalidate();
     (*__iFE)->computeForces(m_force_index);
     );


  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_otherForceTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;

  MSG_DEBUG("Controller::run", "maxCutoff = " << M_SIMULATION->maxCutoff);
  
  // Controller::time_for_parallel = 0;
  // Controller::time_for_parallel1 = 0;
  // Controller::time_for_parallel2 = 0;
  /* --- MAIN LOOP -------------------------------------------------------------------------- */

  for(int current = 0; current < m_timesteps; current++) {
    
    timestep = current*m_dt;
    m_t = timestep;
    
    /* Print out status! */
    if(current%m_statusEvery == 0)
      cout << "@ t = " << timestep << " (" << (double)current*100/m_timesteps << "%), n = "
	   << M_SIMULATION->phase()->returnNofPart() << endl;
    
    /* Print memory usage! */
    /*
      cout << "! memory usage statistics (approximately)." << endl;
      
      cout << "    - cells        "
      << M_MANAGER->cells().size() * (sizeof(Cell) + (11*NUM_NEIGHBORS+2)*sizeof(void*) + sizeof(int))/(1024*1023)
      << " MB" << endl;
      cout << "    - cell links   "
      << M_MANAGER->links().size() * (sizeof(CellLink) + sizeof(void*))/(1024*1024)
      << " MB" << endl;
      cout << "    - particles    "
      << M_SIMULATION->phase()->returnNofPart()*(sizeof(Particle) + 3*sizeof(void*))/(1024*1024)
      << " MB" << endl;
      cout << "    - pairs        " << endl;
      
      for (int ii = 0; ii < M_MANAGER->nColours(); ++ii)
      for (int jj = ii; jj < M_MANAGER->nColours(); ++jj) {
      cout << "          (" << M_MANAGER->species(ii) << ", "
      << M_MANAGER->species(jj) << ")" << endl;
      cout << "             free  "
      << M_MANAGER->cp(ii, jj)->freePairs(0).capacity() * sizeof(Pairdist)/(1024*1024)
      << " MB" << endl;
      cout << "           frozen  "
      << M_MANAGER->cp(ii, jj)->frozenPairs(0).capacity() * sizeof(Pairdist)/(1024*1024)
      << " MB" << endl;
      }
    */

    runSymbols_0();

    // Meters decide themselves if they really have to measure at this timestep
    for(vector<Meter*>::iterator metersIter = M_SIMULATION -> meters() -> begin();
	metersIter != M_SIMULATION -> meters() -> end(); metersIter++) {
      (*metersIter)->measure(timestep);
    }

    integrate();

    //      particle tracking
    //     FOR_EACH_FREE_PARTICLE_C
    //         (M_SIMULATION->phase(), 0,
    //          if(i->mySlot == 62)
    //          {
    //            MSG_DEBUG("Controller::run", i->mySlot << " AFTER integrate(N=" << M_SIMULATION->phase()->returnNofPart() << "):"
    //                << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1] << endl/* << "ie=" << i->tag.doubleByOffset(0)*/);
    //          }
    //         );    
    
    // callables are, e.g., thermostats; they decide themselves, whether activated or not
    FOR_EACH
      (vector<Callable*>,
       (*M_SIMULATION->callables()),
       (*__iFE)->call(current);
       );
    
    /*    computation of derived quantities: currently (2005/04/22) only necessary for 
	  computation of temperature in IntegratorEnergy*/
    FOR_EACH
      (
       list<Node*>, m_children,
       ((Integrator*) *__iFE)->deriveQuantities();
       );
    
  }
  
  /* --- END MAIN LOOP ---------------------------------------------------------------------- */
  
  // measuring of final state
  timestep = m_timesteps*m_dt;
  
  /* Print out status! */
  cout << "@ t = " << timestep << " (100%), n = "
       << M_SIMULATION->phase()->returnNofPart() << endl;

  // additional call of Symbol Calculators
  runSymbols_0();

  for (vector<Meter*>::iterator metersIter = M_SIMULATION -> meters() -> begin();
       metersIter != M_SIMULATION -> meters() -> end(); metersIter++) {
    (*metersIter)->measureNow(timestep);
  }

  MSG_DEBUG("Controller::run", "times needed for\n -pair creation: " << m_pairCreateTime << "\n -pair forces: " << m_pairForceTime << "\n -other forces: " << m_otherForceTime << "\n -time integration: " << m_integrateTime << "\n -initialisation: " << m_initTime << "\n -particle-Symbols: " << m_partSymbolsTime << "\n -bonded Symbols: " << m_bondedSymbolsTime << "\n -non bonded Symbols: " << m_nonBondedSymbolsTime 
// << "\n - special time: " << m_SpecialTime
);

}


Node *Controller::instantiateChild(const string &name)
{
  return Integrator_Factory::byName(name).instantiate(this);
}


void Controller::integrate()
{

  clock_t tInitStart;
  clock_t tInitEnd;
//   time_t tInitStart;
//   time_t tInitEnd;

  tInitStart = clock();
//   std::time(&tInitStart);

  int other_force_index;
  
  Phase* phase = M_SIMULATION->phase();

  /* First integration step */
  FOR_EACH
    (list<Node*>,
     m_children,
     ((Integrator*) *__iFE)->integrateStep1();

    );

//  particle tracking
//     FOR_EACH_FREE_PARTICLE_C
//         (M_SIMULATION->phase(), 0,
//          if(i->mySlot == 180)
//          {
//            MSG_DEBUG("Controller::integrate", i->mySlot << " Nach integrateStep1   N=" << M_SIMULATION->phase()->returnNofPart() << "):"
//                << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1] << endl/* << "ie=" << i->tag.doubleByOffset(0)*/);
//          }
//         );


  /* Recalculate forces. */
  other_force_index = (m_force_index+1)&(FORCE_HIST_SIZE-1);

  // next does NOT clear tag-information, only forces (for velocities) ...
  FOR_EACH_FREE_PARTICLE
    (phase,
     __iSLFE->clear(other_force_index);
     );

  /* Integrators, which use forces in the particle tag, must protect the current forces because
  otherwise they would be deleted by setForNewIntegration() below. At the same time, the
  previous forces are unprotected, because they must be overwritten. Integrators that don't use
  such forces do nothing.*/
  FOR_EACH
      (list<Node*>,
       m_children,
       ((Integrator*) *__iFE)->unprotect(other_force_index);
      );


  // the name is wrong. This clears the unprotected tag, nothing else
  phase->setForNewIntegration();

  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_integrateTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_integrateTime += difftime(tInitEnd, tInitStart);

  tInitStart = clock();
//   std::time(&tInitStart);
      
  FOR_EACH_COLOUR_PAIR
    ( 
     phase->manager(),
     cp->updateConnectedDistances();    
     );

// time_t t0, t1;
// std::time(&t0);
  M_PAIRCREATOR->createDistances();

//   FOR_EACH_COLOUR_PAIR
//     ( 
//      phase->manager(),
//      MSG_DEBUG("Controller::run", "CP(" << cp->firstColour() << cp->secondColour() << "): freePairs=" << cp->freePairs()[0].size() << ", frozenPairs=" << cp->frozenPairs()[0].size());
//      );

  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_pairCreateTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_pairCreateTime += difftime(tInitEnd, tInitStart);

  runSymbols();
  
  // std::time(&t1);
  // MSG_INFO("Controller::run", "2 Parallel/non parallel region took: " << difftime(t1, t0) << " seconds");
//   clock_t start, finish;
//   start = clock();
  
#ifdef _OPENMP 
  double start_time = omp_get_wtime();
#endif
  
  tInitStart = clock();
//   std::time(&tInitStart);

#pragma omp parallel for
  for (int t = 0; t < global::n_threads; ++t) {
    
    vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();
    for(vector<ColourPair*>::iterator __cp = M_MANAGER->colourPairs().begin(); __cp != __end; ++__cp) {
      ColourPair *cp = *__cp;
      vector<GenF*>::iterator pairForcesEnd = cp->pairForces()->end();
      if(cp->freePairsRandom(t).size())
	{
	  // #pragma omp parallel
	  // {
	  
	  // size_t currentPos = 0;
	  // PrimitiveSLEntry<size_t> *i = cp->freePairsRandom(0).first();
	  for (PrimitiveSLEntry<size_t> *i = cp->freePairsRandom(t).first(); i != NULL; i = i->next) {
	    Pairdist* pair = &(cp->freePairs()[t][i->m_val]);	  
	    for (vector<GenF*>::iterator force = cp->pairForces()->begin(); force != pairForcesEnd; ++force) {
	      (*force)->invalidate();
#ifndef _OPENMP
	      (*force)->computeForces(pair, other_force_index);
#else
	      (*force)->computeForces(pair, other_force_index, t);
#endif
	    }
	    
	    // i = i->next;
	    // ++currentPos;
	    //       }
	  }
	  // }
	}
      else {
// 	MSG_DEBUG("Controller::run", "Looping over non-random pair list");
	for (Pairdist *i = cp->freePairs()[t].first(); i != NULL; i = i->next) {
	  
	  Pairdist* pair = i;

//  	  MSG_DEBUG("Controller::run", "pair (" << pair->firstPart()->mySlot << "," << pair->secondPart()->mySlot << ")");
	  
	  for (vector<GenF*>::iterator force = cp->pairForces()->begin(); force != pairForcesEnd; ++force) {
	    (*force)->invalidate();
	    
#ifndef _OPENMP
	    (*force)->computeForces(pair, other_force_index);
#else
	    (*force)->computeForces(pair, other_force_index, t);
#endif
	    
          }
	}
      }
      
      if(cp->frozenPairsRandom(t).size()) {
	for (PrimitiveSLEntry<size_t> *i = cp->frozenPairsRandom(t).first(); i != NULL; i = i->next) {
	  Pairdist* pair = &(cp->frozenPairs()[t][i->m_val]);
	  
	  for (vector<GenF*>::iterator force = cp->pairForces()->begin(); force != pairForcesEnd; ++force) {
	    (*force)->invalidate();
#ifndef _OPENMP
            (*force)->computeForces(pair, other_force_index);
#else
            (*force)->computeForces(pair, other_force_index, t);
#endif
          }
	}
	
      }
      else {
        for (Pairdist *i = cp->frozenPairs()[t].first(); i != NULL; i = i->next) {
          Pairdist* pair = i;
	  
          for (vector<GenF*>::iterator force = cp->pairForces()->begin(); force != pairForcesEnd; ++force) {
            (*force)->invalidate();
#ifndef _OPENMP
	    (*force)->computeForces(pair, other_force_index);
#else
	    (*force)->computeForces(pair, other_force_index, t);
#endif
          }
        }
      }
      
    } // end loop over ColourPairs
  }

#ifdef _OPENMP
  for (int t = 0; t < global::n_threads; t++) {
    list<Node*>::iterator i_begin = m_children.begin();
    list<Node*>::iterator i_end = m_children.end();
    for(list<Node*>::iterator integr = i_begin; integr != i_end; ++integr) {
      size_t _c = ((Integrator*)(*integr))->colour();
      
      FOR_EACH_PARTICLE_C
	(M_SIMULATION->phase(), _c,
	 // MSG_DEBUG("Controller::run", " integr = " << ((Integrator*)(*integr))->className());
	 ((Integrator*)(*integr))->mergeCopies(__iSLFE, t, other_force_index);
	 );
    }
  }
#endif
  
#ifdef _OPENMP
  double end_time = omp_get_wtime(); 
  double wtick = omp_get_wtick();
  
  // MSG_DEBUG("Controller::run", " time elapsed = " << end_time - start_time);
  // MSG_DEBUG("Controller::run", " ticks elapsed = " << wtick);
#endif
  
  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_pairForceTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_pairForceTime += difftime(tInitEnd, tInitStart);

  tInitStart = clock();
//   std::time(&tInitStart);
  
  for(size_t c = 0; c < M_PHASE->nColours(); ++c) {   
    for (Particle *i = M_PHASE->particles(c).first(); i != NULL; i = i->next) {    
      
      for (vector<GenF*>::iterator force = (*M_SIMULATION->particleForces())[c]->begin(); force != (*M_SIMULATION->particleForces())[c]->end(); ++force) {
	(*force)->invalidate();
	(*force)->computeForces(i, other_force_index);
      }
      
    }
  }
  
  
  FOR_EACH
    (vector<GenF*>,
     (*M_SIMULATION->otherForces()),
     (*__iFE)->invalidate();
     (*__iFE)->computeForces(other_force_index);
     );
  
  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_otherForceTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_otherForceTime += difftime(tInitEnd, tInitStart);

  tInitStart = clock();
//   std::time(&tInitStart);

  (++m_force_index) &= (FORCE_HIST_SIZE-1);

  /* Second integration step */
  FOR_EACH
    (list<Node*>,
     m_children,
     ((Integrator*) *__iFE)->integrateStep2();
    );

  tInitEnd = clock();
//   std::time(&tInitEnd);
  m_integrateTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//   m_integrateTime += difftime(tInitEnd, tInitStart);

}


void Controller::runSymbols() {
//    MSG_DEBUG("Controller::runSymbols", "START");

  clock_t tInitStart;
  clock_t tInitEnd;
//   time_t tInitStart;
//   time_t tInitEnd;

  Phase* phase = M_PHASE;
  size_t stage = 0;
  bool cpsFinished = false;

  while(stage <= Particle::s_maxStage/*PCA_MAX_STAGE*/ || !cpsFinished)
    {
//       MSG_DEBUG("Controller::runSymbols", "LOOP-START: stage = " << stage << ", Pstage = " << Particle::s_maxStage << ", cpsFinished = " << cpsFinished);
      /* Run ParticleCaches... Means calculator information that 
	 is particle dependent. Can rely on calculator outputs from
	 same and lower stages and on cache outputs from lower stages */
      //      if(stage <= Particle::s_maxStage/*PCA_MAX_STAGE*/)
      //      {
//       ManagerCell* manager = phase->manager();

      tInitStart = clock();
//       std::time(&tInitStart);

      size_t nCols = phase->nColours();
      for (size_t col = 0; col < nCols; ++col) 
	{
	  if (/* !Particle::s_cached_properties[col][stage].empty()*/
	      Particle::s_cached_properties[col].size() > stage) 
	    {
	      
	      vector<ParticleCache*>::iterator __begin, __end;
	      
	      __begin = Particle::s_cached_properties[col][stage].begin();
	      __end = Particle::s_cached_properties[col][stage].end();
	      
	      FOR_EACH_FREE_PARTICLE_C
		(
		 phase,
		 col,
		 for (vector<ParticleCache*>::iterator pc = __begin; pc != __end; ++pc)
		   (*pc)->computeCacheFor(__iSLFE);
		 );
	    }
	}
      //       }

      tInitEnd = clock();
//       std::time(&tInitEnd);
      m_partSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//       m_partSymbolsTime += difftime(tInitEnd,tInitStart);

      tInitStart = clock();
//       std::time(&tInitStart);
            
      // the bonded triplets
      //FIXME: parallelise!
      size_t tlSize = phase->tripletLists()->size(); 
      for(size_t itl = 0; itl < tlSize; ++itl) {
	if(phase->maxBondedStage(itl) >= stage) {
	  tripletList* tL = phase->returnTripletList(itl);
	  // loop over triplets of current connected list 
	  tripletListItr trEnd = tL->end();
	  for(tripletListItr tr = tL->begin(); tr != trEnd; ++tr) {
	    tr->runBondedTripletCalculators(stage, itl);
	  }
	}
      }

      tInitEnd = clock();
//       std::time(&tInitEnd);
      m_bondedSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//       m_bondedSymbolsTime += difftime(tInitEnd,tInitStart);


#ifndef _OPENMP
// time(&t0);
// t0 = clock();
      cpsFinished = true;
      FOR_EACH_COLOUR_PAIR
	(M_MANAGER,
	 // Run all calculators for bonded pairs in current stage.
	 // FIXME: parallelise!
	 // loop over vector of connected lists

	 tInitStart = clock();
// 	 std::time(&tInitStart);

	 size_t clSize = cp->connectedLists()->size(); 
	 for(size_t icl = 0; icl < clSize; ++icl) {
	   if(cp->maxBondedStage(icl) >= stage) {
	     PairList* pL = cp->connectedList(icl);
	     // loop over pairs of current connected list 
	     for(Pairdist *pair = pL->first(); 
		 pair != NULL; pair = pair->next) {
	       pair->runBondedPairCalculators(stage, icl);
	     }
	     if (cp->maxBondedStage(icl) == stage) cpsFinished = cpsFinished && true; 
	     else cpsFinished = false;
	   }
	 }
//       MSG_DEBUG("Controller::runSymbols", "IN CP: before non-bonded");

	 tInitEnd = clock();
// 	 std::time(&tInitEnd);
	 m_bondedSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
// 	 m_bondedSymbolsTime += difftime(tInitEnd,tInitStart);

	 tInitStart = clock();
// 	 std::time(&tInitStart);

// 	 clock_t tSpecial;
// 	 clock_t tSpecial2; 
// 	 tSpecial = clock();

	 /* Run all calculators for non-bonded pairs in current stage. */
	 if (cp->maxStage() >= stage) {
//  	 MSG_DEBUG("Controller::run", "nonBonded-calc: stage-check TRUE; stage=" << stage);
	   FOR_EACH
	     (vector<ValCalculator*>,
	      cp->valCalculators(stage),
// 	      MSG_DEBUG("Controller::run", "nonBonded-calc: " << (*__iFE)->className() << ", stage=" << stage);
	      FOR_EACH_PAIR__PARALLEL
	      (Controller,
	       cp,
 	       (*__iFE)->compute(pair);
 	       );
	      // old version (2013-05-22); if new version works, check 
	      // if you can remove this member function of Pairdist
	      // 	      pair->runCalculatorsForStage(stage);
	      );
	   
	   if (cp->maxStage() == stage) cpsFinished = cpsFinished && true; 
	   else cpsFinished = false;
	 }

// 	 tSpecial2 = clock();
// 	 m_SpecialTime += (tSpecial2 - tSpecial)/(double)CLOCKS_PER_SEC;
	 
	 tInitEnd = clock();
// 	 std::time(&tInitEnd);
	 m_nonBondedSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
// 	 m_nonBondedSymbolsTime += difftime(tInitEnd,tInitStart);

	 );
      ++stage;

#else

     cpsFinished = true;
     
     // Run all calculators for bonded pairs in current stage.
     // FIXME: parallelise! Currently NOT parallel!!!     
     FOR_EACH_COLOUR_PAIR
       (M_MANAGER,
	// next line commented out because maxStage() and maxBondedStage decoupled now (2010-01-18)
// 	if (cp->maxStage() >= stage) {

	  // loop over vector of connected lists
	  size_t clSize = cp->connectedLists()->size(); 
	  for(size_t icl = 0; icl < clSize; ++icl) {
	    if(cp->maxBondedStage(icl) >= stage) {
	      PairList* pL = cp->connectedList(icl);
	      // loop over pairs of current connected list 
	      for(Pairdist *pair = pL->first(); 
		  pair != NULL; pair = pair->next) {
		pair->runBondedPairCalculators(stage, icl, 0/*thread number when parallelised!!!*/);
	      }
	      if (cp->maxBondedStage(icl) == stage) cpsFinished = cpsFinished && true; 
	      else cpsFinished = false;
	    }
	  }
// 	} // end of if (cp->maxStage() >= stage)
	);


// time(&t0);
// t0 = clock();
// vector<ColourPair*>::iterator cp;
//           vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
//           for(vector<ColourPair*>::iterator cp = M_MANAGER->colourPairs().begin(); cp != __end; ++cp) {                                                
//             if ((*cp)->maxStage() >= stage) {
#pragma omp parallel for
        for (int t = 0; t < global::n_threads; ++t) {
          vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
          for(vector<ColourPair*>::iterator cp = M_MANAGER->colourPairs().begin(); cp != __end; ++cp) {                                            
//             ColourPair *cp = *__cp; 
        
            if ((*cp)->maxStage() >= stage) {
   
	      // Run all calculators for non-bonded pairs        
              if((*cp)->freePairsRandom(t).size())
              {
                for (PrimitiveSLEntry<size_t> *i = (*cp)->freePairsRandom(t).first(); i != NULL; i = i->next) {
                  Pairdist* pair = &((*cp)->freePairs()[t][i->m_val]);
		  // FIXME: it should be possible to directly write the contents of the called
		  // function Pairdist::runCalculatorsForStage(size_t stage, int thread_no). Check
		  // the same for other similar places
                  pair->runCalculatorsForStage(stage, t);
                }
              }  
              else 
		{
		  for (Pairdist *i = (*cp)->freePairs()[t].first(); i != NULL; i = i->next) {
		    Pairdist* pair = i;
		    pair->runCalculatorsForStage(stage, t);
		  }
		}
	      
              if((*cp)->frozenPairsRandom(t).size()) 
		{
		  for (PrimitiveSLEntry<size_t> *i = (*cp)->frozenPairsRandom(t).first(); i != NULL; i = i->next) {
		    Pairdist* pair = &((*cp)->frozenPairs()[t][i->m_val]);
		    
		    pair->runCalculatorsForStage(stage, t);
		  }
		}
              else 
		{
		  for (Pairdist *i = (*cp)->frozenPairs()[t].first(); i != NULL; i = i->next) {
		    Pairdist* pair = i;
		    
		    pair->runCalculatorsForStage(stage, t);
		  }
		}
	      
              if ((*cp)->maxStage() == stage) cpsFinished = cpsFinished && true;
              else cpsFinished = false;
            }
          }
        }
// time(&t2); 
// t2 = clock();
// Add the copies alltogether to get the originally calculated value, saved in the original Particle slot
        for (size_t t = 0; t < global::n_threads; ++t) {
          vector<ColourPair*>::iterator __begin = M_MANAGER->colourPairs().begin();
          vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
          for(vector<ColourPair*>::iterator cp = __begin; cp != __end; ++cp) {

            if ((*cp)->maxStage() >= stage) {

              vector<ValCalculator*>::iterator vc_begin = (*cp)->valCalculatorParts(stage).begin();
              vector<ValCalculator*>::iterator vc_end = (*cp)->valCalculatorParts(stage).end();

              while (vc_begin != vc_end) {
                
                ((ValCalculatorPart*)(*vc_begin))->mergeCopies((*cp), t);

                ++ vc_begin;
              }
            }
          }
        }

        ++stage;
// time(&t1);
// t1 = clock();

#endif

// Controller::time_for_parallel += /*double(t1-t0)/CLOCKS_PER_SEC*/difftime(t1, t0);
// Controller::time_for_parallel1 += /*double(t2-t0)/CLOCKS_PER_SEC*/difftime(t2, t0);
// Controller::time_for_parallel2 += /*double(t2-t1)/CLOCKS_PER_SEC*/difftime(t2, t1);
//     MSG_INFO("COntroller::run", "Parallel/non parallel region took: " << /*double(*/t1-t0/*)/CLOCKS_PER_SEC*/ << " seconds. all sim parallel = " << Controller::time_for_parallel << " all sim only parallel part = " << Controller::time_for_parallel1 /* << " only parallel part = " << int(t2-t0)/CLOCKS_PER_SEC << " adding the particles together = " << int(t2-t1)/CLOCKS_PER_SEC*/);


    }
//    MSG_DEBUG("Controller::runSymbols", "END");
}

void Controller::runSymbols_0() {
//   MSG_DEBUG("Controller::runSymbols_0", "START");

  clock_t tInitStart;
  clock_t tInitEnd;
//   time_t tInitStart;
//   time_t tInitEnd;

  Phase* phase = M_PHASE;
  size_t stage = 0;
  bool cpsFinished = false;
  while(stage <= Particle::s_maxStage_0/*PCA_MAX_STAGE*/ || !cpsFinished)
    {
    //MSG_DEBUG("Controller::runSymbols_0", "stage = " << Particle::s_maxStage_0);
    
      /* Run ParticleCaches... Means calculator information that 
	 is particle dependent. Can rely on calculator outputs from
	 same and lower stages and on cache outputs from lower stages */
      //      if(stage <= Particle::s_maxStage_0/*PCA_MAX_STAGE*/)
      //      {
//       ManagerCell* manager = phase->manager();

      tInitStart = clock();
//       std::time(&tInitStart);

      size_t nCols = phase->nColours();
      for (size_t col = 0; col < nCols; ++col) 
	{
	  if (/* !Particle::s_cached_properties_0[col][stage].empty()*/
	      Particle::s_cached_properties_0[col].size() > stage) 
	    {
	      
	      vector<ParticleCache*>::iterator __begin, __end;
	      
	      __begin = Particle::s_cached_properties_0[col][stage].begin();
	      __end = Particle::s_cached_properties_0[col][stage].end();
	      
	      FOR_EACH_FREE_PARTICLE_C
		(
		 phase,
		 col,
		 for (vector<ParticleCache*>::iterator pc = __begin; pc != __end; ++pc)
		   (*pc)->computeCacheFor(__iSLFE);
		 );
	    }
	}
      //       }

      tInitEnd = clock();
//       std::time(&tInitEnd);
      m_partSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//       m_partSymbolsTime += difftime(tInitEnd, tInitStart);

      tInitStart = clock();
//       std::time(&tInitStart);
      
      // the bonded triplets
      //FIXME: parallelise!
      size_t tlSize = phase->tripletLists()->size(); 
 //    MSG_DEBUG("Controller::runSymbols_0", "tlSize = " << tlSize); 
      for(size_t itl = 0; itl < tlSize; ++itl) {
	    if(phase->maxBondedStage_0(itl) >= stage) {
	      tripletList* tL = phase->returnTripletList(itl);
	      // loop over pairs of current connected list 
	      tripletListItr trEnd = tL->end();
	      
	      for(tripletListItr tr = tL->begin(); tr != trEnd; ++tr) {
	        tr->runBondedTripletCalculators_0(stage, itl);
	      }
	    }
      }


      tInitEnd = clock();
//       std::time(&tInitEnd);
      m_bondedSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//       m_bondedSymbolsTime += difftime(tInitEnd, tInitStart);

      /* Run all calculators in current stage. These are calculators
	 like the shear rate calculator which depend on information
	 from lower stages.
      */

#ifndef _OPENMP
// time(&t0);
// t0 = clock();
      cpsFinished = true;
      FOR_EACH_COLOUR_PAIR
	(M_MANAGER,
	 // Run all calculators for bonded pairs in current stage.
	 // FIXME: parallelise!
	 // loop over vector of connected lists

	 tInitStart = clock();
// 	 std::time(&tInitStart);

	 size_t clSize = cp->connectedLists()->size(); 
	 for(size_t icl = 0; icl < clSize; ++icl) {
	   if(cp->maxBondedStage_0(icl) >= stage) {
	     PairList* pL = cp->connectedList(icl);
	     // loop over pairs of current connected list 
	     for(Pairdist *pair = pL->first(); 
		 pair != NULL; pair = pair->next) {
	       pair->runBondedPairCalculators_0(stage, icl);
	     }
	     if (cp->maxBondedStage_0(icl) == stage) cpsFinished = cpsFinished && true; 
	     else cpsFinished = false;
	   }
	 }

      tInitEnd = clock();
//       std::time(&tInitEnd);
      m_bondedSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//       m_bondedSymbolsTime += difftime(tInitEnd, tInitStart);

	 tInitStart = clock();
// 	 std::time(&tInitStart);
	 
	 /* Run all calculators for non-bonded pairs in current stage. */
	 if (cp->maxStage_0() >= stage) {
	   FOR_EACH_PAIR__PARALLEL
	     (Controller,
	      cp,
	      // new version (2013-05-22)
 	      FOR_EACH
 	      (vector<ValCalculator*>,
 	       cp->valCalculators_0(stage),
 	       (*__iFE)->compute(pair);
 	       );
	      // old version (2013-05-22); if new version works, check 
	      // if you can remove this member function of Pairdist
// 	      pair->runCalculatorsForStage_0(stage);
	);

      tInitEnd = clock();
//       std::time(&tInitEnd);
      m_nonBondedSymbolsTime += (tInitEnd - tInitStart)/(double)CLOCKS_PER_SEC;
//       m_nonBondedSymbolsTime += difftime(tInitEnd, tInitStart);

	   if (cp->maxStage_0() == stage) cpsFinished = cpsFinished && true; 
	   else cpsFinished = false;
	 }
	 );
      ++stage;
// time(&t1);
// t1 = clock();


#else // parallel case
      cpsFinished = true;

     // Run all calculators for bonded pairs in current stage.
     // FIXME: parallelise! Currently NOT parallel!!!     
     FOR_EACH_COLOUR_PAIR
       (M_MANAGER,
	// next line commented out because now (2010-01-18) maxStage_0() and maxBondedStage_0() have been decoupled
// 	if (cp->maxStage_0() >= stage) {
	  // loop over vector of connected lists
	  size_t clSize = cp->connectedLists()->size(); 
//	  MSG_DEBUG("Controller::runSymbols_0", "clSize = " << clSize); 
	  for(size_t icl = 0; icl < clSize; ++icl) {
	    if(cp->maxBondedStage_0(icl) >= stage) {
	      PairList* pL = cp->connectedList(icl);
	      // loop over pairs of current connected list 
	      for(Pairdist *pair = pL->first(); pair != NULL; pair = pair->next) {
		    pair->runBondedPairCalculators_0(stage, icl, 0/*thread number when parallelised!!!*/);
	      }
	      if (cp->maxBondedStage_0(icl) == stage) cpsFinished = cpsFinished && true; 
	      else cpsFinished = false;
	    }
	  }
	//} // end of if (cp->maxStage_0() >= stage)
	);




// vector<ColourPair*>::iterator cp;
//           vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
//           for(vector<ColourPair*>::iterator cp = M_MANAGER->colourPairs().begin(); cp != __end; ++cp) {                                                
//             if ((*cp)->maxStage_0() >= stage) {
#pragma omp parallel for 
        for (int t = 0; t < global::n_threads; ++t) {
          vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
          for(vector<ColourPair*>::iterator cp = M_MANAGER->colourPairs().begin(); cp != __end; ++cp) {                                            
//             ColourPair *cp = *__cp; 
            if ((*cp)->maxStage_0() >= stage) {
              if((*cp)->freePairsRandom(t).size())
              {
                for (PrimitiveSLEntry<size_t> *i = (*cp)->freePairsRandom(t).first(); i != NULL; i = i->next) {
                  Pairdist* pair = &((*cp)->freePairs()[t][i->m_val]);
            
                  pair->runCalculatorsForStage_0(stage, t);
                }
              }  
              else 
              {
                for (Pairdist *i = (*cp)->freePairs()[t].first(); i != NULL; i = i->next) {
                  Pairdist* pair = i;
                  pair->runCalculatorsForStage_0(stage, t);
                }
              }
          
              if((*cp)->frozenPairsRandom(t).size()) 
              {
                for (PrimitiveSLEntry<size_t> *i = (*cp)->frozenPairsRandom(t).first(); i != NULL; i = i->next) {
                  Pairdist* pair = &((*cp)->frozenPairs()[t][i->m_val]);
            
                  pair->runCalculatorsForStage_0(stage, t);
                }
              }
              else 
              {
                for (Pairdist *i = (*cp)->frozenPairs()[t].first(); i != NULL; i = i->next) {
                  Pairdist* pair = i;
          
                  pair->runCalculatorsForStage_0(stage, t);
                }
              }
              if ((*cp)->maxStage_0() == stage) cpsFinished = cpsFinished && true;
              else cpsFinished = false;
            }
          }
        }
// time(&t2); 
// t2 = clock();
// Add the copies alltogether to get the originally calculated value, saved in the original Particle slot
        for (size_t t = 0; t < global::n_threads; ++t) {
          vector<ColourPair*>::iterator __begin = M_MANAGER->colourPairs().begin();
          vector<ColourPair*>::iterator __end = M_MANAGER->colourPairs().end();   
          for(vector<ColourPair*>::iterator cp = __begin; cp != __end; ++cp) {

            if ((*cp)->maxStage_0() >= stage) {

              vector<ValCalculator*>::iterator vc_begin = (*cp)->valCalculatorParts(stage).begin();
              vector<ValCalculator*>::iterator vc_end = (*cp)->valCalculatorParts(stage).end();

              while (vc_begin != vc_end) {
                
                ((ValCalculatorPart*)(*vc_begin))->mergeCopies((*cp), t);

                ++ vc_begin;
              }
            }
          }
        }

        ++stage;
// time(&t1);
// t1 = clock();
#endif
// Controller::time_for_parallel += /*double(t1-t0)/CLOCKS_PER_SEC*/difftime(t1, t0);
// Controller::time_for_parallel1 += /*double(t2-t0)/CLOCKS_PER_SEC*/difftime(t2, t0);
// Controller::time_for_parallel2 += /*double(t2-t1)/CLOCKS_PER_SEC*/difftime(t2, t1);
//     MSG_INFO("LinkedListCreator::createDistances", "Parallel/non parallel region took: " << /*double(*/t1-t0/*)/CLOCKS_PER_SEC*/ << " seconds. all sim parallel = " << Controller::time_for_parallel << " all sim only parallel part = " << Controller::time_for_parallel1 /* << " only parallel part = " << int(t2-t0)/CLOCKS_PER_SEC << " adding the particles together = " << int(t2-t1)/CLOCKS_PER_SEC*/);


  }
}


void Controller::setup() {
  NodeManyChildren::setup();

//   MSG_DEBUG("Controller::setup", "Valgrindcheck *this: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(*this));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_timesteps: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_timesteps));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_dt: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_dt));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_t: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_t));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_statusEvery: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_statusEvery));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_force_index: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_force_index));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_toSetupAfterParticleCreation: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_toSetupAfterParticleCreation));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck m_toSetupAfterParticleCreation.begin: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(m_toSetupAfterParticleCreation.begin()));
//   MSG_DEBUG("Controller::setup", "Valgrindcheck (*(m_toSetupAfterParticleCreation.begin())): "<< VALGRIND_CHECK_VALUE_IS_DEFINED((*(m_toSetupAfterParticleCreation.begin()))));

}


void Controller::registerForSetupAfterParticleCreation (Node* callable)
{

  m_toSetupAfterParticleCreation.push_back(callable); 

}
