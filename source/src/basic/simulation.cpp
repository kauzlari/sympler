/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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



#include <fstream>
#include <algorithm>

#include "simulation.h"
#include "manager_cell.h"
#include "wall_triangle.h"
#include "reflector_stochastic.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "pair_creator.h"
#include "val_calculator_part.h"

#ifdef _OPENMP
  #include "omp.h"
#endif

// #include "valgrind/memcheck.h"


using namespace std;

size_t global::n_threads;
// vector<int> Simulation::nDoublesProCol;


//---- Constructors/Destructor ----

Simulation::Simulation()
  : NodeManyChildren(), m_input_from_results(false), m_controller(NULL), m_phase(NULL)
{
  init();
}


Simulation::~Simulation()
{
  /* m_phase, meters and forces will be deleted by the destructor of
     NodeManyChildren */
}



//---- Methods ----

void Simulation::setSymbolStages()
{
  bool finished = false;
  size_t counter = 1;
  while(!finished)
  {
    MSG_DEBUG("Simulation::setSymbolStages", "iteration #" << counter);
    if(counter > m_stageSteps)
      throw gError("Simulation::setSymbolStages", "I was not able to determine all Symbol stages in " + ObjToString(counter) + " iterations. Check for inconsistencies in your runtime-compiled expressions or increase 'stageIterations'.");
    finished = true;
    finished = phase()->findStages() && finished;
    finished = phase()->findStages_0() && finished;
    FOR_EACH_COLOUR_PAIR
    (
      phase()->manager(),
      finished = cp->findStages() && finished;
      finished = cp->findStages_0() && finished;
    );

      // !!! findStage() part must come first, because, if first argument is false, the second part will not be evaluated at all !!!
    for(vector<ParticleCache*>::iterator i = Particle::s_cached_flat_properties.begin(); i != Particle::s_cached_flat_properties.end(); ++i)
    {
      finished = (*i)->findStage() && finished;
    }
    for(vector<ParticleCache*>::iterator i = Particle::s_cached_flat_properties_0.begin(); i != Particle::s_cached_flat_properties_0.end(); ++i)
    {
      finished = (*i)->findStage_0() && finished;
    }
    ++counter;
  }
}

void Simulation::sortSymbolsByStages() {
  phase()->sortStages();
  phase()->sortStages_0();
  FOR_EACH_COLOUR_PAIR
  (
   phase()->manager(),
   cp->sortStages();
   cp->sortStages_0();
   );
  
  Particle::sortStages();
  Particle::sortStages_0();
}

#ifdef _OPENMP
// Set the offset to the copy vectors for the parallel version
// Set the proper places to write in the copy vectors for each ValCalculator
void Simulation::setupCopyVectors()
{
  size_t stage = 0;
  bool cpsFinished = false;
  while(!cpsFinished) {
    /*      size_t temp_cp = 0;*/
    vector<size_t> temp_cp(m_phase->manager()->nColours(),0);
    cpsFinished = true;
    FOR_EACH_COLOUR_PAIR
      (phase()->manager(),
       if (cp->maxStage() >= stage){

	 for(vector<ValCalculator*>::iterator i = cp->valCalculatorParts(stage).begin(); i != cp->valCalculatorParts(stage).end(); ++i)
	   {
	     ((ValCalculatorPart*)(*i))->vectorSlots().first = temp_cp[cp->firstColour()];
	     ((ValCalculatorPart*)(*i))->vectorSlots().second = temp_cp[cp->secondColour()];

                temp_cp[cp->firstColour()] += ((ValCalculatorPart*)(*i))->setNumOfDoubles();
		// do the next only if different colours
		if(cp->firstColour() != cp->secondColour()) 
		  temp_cp[cp->secondColour()] += ((ValCalculatorPart*)(*i))->setNumOfDoubles();

                string tempStrFirst = "copy" + ObjToString(cp->firstColour());
                string tempStrSecond = "copy" + ObjToString(cp->secondColour());

		// The vector doubles have to be filled with 0s at the end of every stage

                for (size_t t = 0; t < global::n_threads; t++) {
                  if (!Particle::s_tag_format[cp->firstColour()].attrExists(tempStrFirst + ObjToString(t))) {
		            throw gError("Simulation::setupCopyVectors", "copy-vector \"" + tempStrFirst + ObjToString(t) + "\" not found! Contact a programmer.");
                  }
                    if (Particle::s_tag_format[cp->firstColour()].vectorDoubleByName(tempStrFirst + ObjToString(t))->size() < temp_cp[cp->firstColour()])
                      Particle::s_tag_format[cp->firstColour()].vectorDoubleByName(tempStrFirst + ObjToString(t))->resize(temp_cp[cp->firstColour()]);

		    if (!Particle::s_tag_format[cp->secondColour()].attrExists(tempStrSecond + ObjToString(t))) {
		      throw gError("Simulation::setupCopyVectors", "copy-vector \"" + tempStrSecond + ObjToString(t) + "\" not found! Contact a programmer.");

		    }
		    if (Particle::s_tag_format[cp->secondColour()].vectorDoubleByName(tempStrSecond + ObjToString(t))->size() < temp_cp[cp->secondColour()]) {
                      Particle::s_tag_format[cp->secondColour()].vectorDoubleByName(tempStrSecond + ObjToString(t))->resize(temp_cp[cp->secondColour()]);
		      
		    }
                 
		  ((ValCalculatorPart*)(*i))->copySlots()[t].first = Particle::s_tag_format[cp->firstColour()].attrByName(tempStrFirst + ObjToString(t)).offset; 
		  ((ValCalculatorPart*)(*i))->copySlots()[t].second = Particle::s_tag_format[cp->secondColour()].attrByName(tempStrSecond + ObjToString(t)).offset;
		  
            }      
          }

          // the CPs may have different max stages
          if (cp->maxStage() == stage) cpsFinished = cpsFinished && true;
          else cpsFinished = false;
        }
      );
      ++stage;
    }


  // Set the offset to the copy vectors for the parallel version for the stage 0
  // Set the proper places to write in the copy vectors for each ValCalculator
  
  stage = 0;
  cpsFinished = false;
  while(!cpsFinished)
    {
      vector<size_t> temp_cp(m_phase->manager()->nColours(),0);
      cpsFinished = true;
      // loop over ColourPair's
      FOR_EACH_COLOUR_PAIR
          (phase()->manager(),

        if (cp->maxStage_0() >= stage){

	  // loop over ValCalculator's
          for(vector<ValCalculator*>::iterator i = cp->valCalculatorParts_0(stage).begin(); i != cp->valCalculatorParts_0(stage).end(); ++i)
            {
	      ((ValCalculatorPart*)(*i))->vectorSlots().first = temp_cp[cp->firstColour()];
	      ((ValCalculatorPart*)(*i))->vectorSlots().second = temp_cp[cp->secondColour()];

	      temp_cp[cp->firstColour()] += ((ValCalculatorPart*)(*i))->setNumOfDoubles();
	      temp_cp[cp->secondColour()] += ((ValCalculatorPart*)(*i))->setNumOfDoubles();

	      string tempStrFirst = "copy" + ObjToString(cp->firstColour());
	      string tempStrSecond = "copy" + ObjToString(cp->secondColour());
      
	      // The vector doubles have to be filled with 0s at the end of every stage   
	      // Loop over the number of threads
	      for (size_t t = 0; t < global::n_threads; t++) {
		
		    if (!Particle::s_tag_format[cp->firstColour()].attrExists(tempStrFirst + ObjToString(t))) {
		      throw gError("Simulation::setupCopyVectors", "copy-vector \"" + tempStrFirst + ObjToString(t) + "\" not found! Contact a programmer.");
          }
                    if (Particle::s_tag_format[cp->firstColour()].vectorDoubleByName(tempStrFirst + ObjToString(t))->size() < temp_cp[cp->firstColour()])
                      Particle::s_tag_format[cp->firstColour()].vectorDoubleByName(tempStrFirst + ObjToString(t))->resize(temp_cp[cp->firstColour()]);

                  if (!Particle::s_tag_format[cp->secondColour()].attrExists(tempStrSecond + ObjToString(t))) {
        		    throw gError("Simulation::setupCopyVectors", "copy-vector \"" + tempStrSecond + ObjToString(t) + "\" not found! Contact a programmer.");
                  }
                    if (Particle::s_tag_format[cp->secondColour()].vectorDoubleByName(tempStrSecond + ObjToString(t))->size() < temp_cp[cp->secondColour()])
                      Particle::s_tag_format[cp->secondColour()].vectorDoubleByName(tempStrSecond + ObjToString(t))->resize(temp_cp[cp->secondColour()]);

                    ((ValCalculatorPart*)(*i))->copySlots()[t].first = Particle::s_tag_format[cp->firstColour()].attrByName("copy" + ObjToString(cp->firstColour()) + ObjToString(t)).offset;
                    ((ValCalculatorPart*)(*i))->copySlots()[t].second = Particle::s_tag_format[cp->secondColour()].attrByName("copy" + ObjToString(cp->secondColour()) + ObjToString(t)).offset;

                } // end of loop over threads

            } // end of loop over ValCalculator's

          // the CPs may have different max stages
          if (cp->maxStage_0() == stage) cpsFinished = cpsFinished && true;
          else cpsFinished = false;
        }
      ); // end of loop over ColourPair's
      ++stage;
    }

    // Setting the offsets for the force-copies in Particle tag. Here, the memory space already reserved by the ValCalculators
    // (i.e. tag space) will be reused and extended, if necessary. No need to save brand new space for the Force copies!

    vector<size_t> temp_nd(m_phase->manager()->nColours(),0);

    // loop over the Integrators
    list<Node*>::iterator i_begin = m_controller->integrators()->begin();
    list<Node*>::iterator i_end = m_controller->integrators()->end();
    for(list<Node*>::iterator integr = i_begin; integr != i_end; ++integr) {
	size_t _c = ((Integrator*)(*integr))->colour();
	((Integrator*)(*integr))->posInVec() = temp_nd[_c];

	temp_nd[_c] += ((Integrator*)(*integr))->numCopyDoubles();

	// loop over the threads used
	for (size_t t = 0; t < global::n_threads; t++) {
	  if (Particle::s_tag_format[_c].attrExists("copy" + ObjToString(_c) + ObjToString(t))) {
	    if (Particle::s_tag_format[_c].vectorDoubleByName("copy" + ObjToString(_c) + ObjToString(t))->size() < temp_nd[_c])
	      Particle::s_tag_format[_c].vectorDoubleByName("copy" + ObjToString(_c) + ObjToString(t))->resize(temp_nd[_c]);

	  }
	  else {
	    throw gError("Simulation::setupCopyVectors", "copy-vector \"copy" + ObjToString(_c) + ObjToString(t) + "\" not found! Contact a programmer.");

	  }
	  
	  ((Integrator*)(*integr))->offsetToVec()[t] = Particle::s_tag_format[_c].attrByName("copy" + ObjToString(_c) + ObjToString(t)).offset;

	  vector<ColourPair*>::iterator __end = m_phase->manager()->colourPairs().end();   
	  vector<ColourPair*>::iterator __begin = m_phase->manager()->colourPairs().begin();
	  for(vector<ColourPair*>::iterator __cp = __begin; __cp != __end; ++__cp) { 
	    for (vector<GenF*>::iterator _f = (*__cp)->pairForces()->begin(); _f != (*__cp)->pairForces()->end(); ++_f) {
	      (*_f)->setForceSlots(((Integrator*)(*integr)), t);
	    }
	  }
	  
	} // end of loop over threads
  } // end of loop over the Integrators
       
}
#endif

void Simulation::setupMeters()
{
  for (vector<Meter*>::iterator i = m_meters.begin(); i != m_meters.end(); ++i) {
    (*i)->aboutToStart();
  }
}


void Simulation::run()
{
  setup();

  MSG_DEBUG("Simulation::run", "The following species are registered:");

  for (size_t i = 0; i < m_phase->manager()->nColours(); i++) {
    cout << m_phase->manager()->species(i) << endl;
  }


  MSG_DEBUG("Simulation::run", "The following forces are registered:");

  vector<ColourPair*>::iterator __end = m_phase->manager()->colourPairs().end();
  for(vector<ColourPair*>::iterator __cp = m_phase->manager()->colourPairs().begin(); __cp != __end; ++__cp) {
    for (vector<GenF*>::iterator i = (*__cp)->pairForces()->begin(); i != (*__cp)->pairForces()->end(); i++) {
      cout << (*i)->className() << endl;
      cout << "... which is a pair force" << endl;
    }
  }
  for (size_t col = 0; col < m_phase->manager()->nColours(); ++col) {
    for (vector<GenF*>::iterator i = m_particle_forces[col]->begin(); i != m_particle_forces[col]->end(); i++) {
      cout << (*i)->className() << endl;
    }
  }
  for (vector<GenF*>::iterator i = m_other_forces.begin(); i != m_other_forces.end(); i++) {
    cout << (*i)->className() << endl;
    cout << "... which is an \"other\" force" << endl;
  }

  MSG_INFO("Simulation::run", "Starting Simulation...");
  m_controller->run();
  MSG_INFO("Simulation::run", "End of Simulation.");

  FOR_EACH
    (vector<Meter*>,
     m_meters,
     (*__iFE)->flush();
    );

  if(m_input_from_results)
    writeWhenDone();
}


void Simulation::readWithArg(const int& argc, char* argv[])
{
  xmlDoc *doc = NULL;
  xmlNode *root_element = NULL;

  if(argc != 2)
    throw gError("Simulation::read", "Wrong number of arguments.\nSyntax is\n\n\t./sympler "
                 "INPUT\n\nwhere INPUT is your input file.\n");
/*  ifstream in(argv[1]);
    parser p1(in);
    SGML e1;
    // check if input file exists
    if(!in) {
		cout << "File not found." << endl;
    throw gError("Simulation::read", "Could not read input file '" +
    ObjToString(argv[1]) + "'.");
    }
    p1.GetSGML(e1);
    NodeManyChildren::read(e1);*/

   LIBXML_TEST_VERSION

  /* Parse the file and get the DOM */
  doc = xmlReadFile(argv[1], NULL, 0);

  if (doc == NULL) {
    throw gError("Simulation::read", "Could not parse input file: " + string(argv[1]));
  }

  root_element = xmlDocGetRootElement(doc);

  NodeManyChildren::read(root_element);

  xmlFreeDoc(doc);

  xmlCleanupParser();
}


void Simulation::init()
{
  m_properties.setClassName("Simulation");

  m_properties.setDescription(
			      "This is the base object which sets global simulation parameters "
			      "and contains the Forces, a Phase and the Meters as its children."
			      );

  STRINGPC
    (simName, m_name,
     "Name of the simulation. Used for filename generation.");

#ifdef _OPENMP
  INTPC
    (nThreads, global::n_threads, 0,
     "Number of parallel threads to use in this calculation."
     "Please note that the optimum is: nThreads </= nProcessors.");
#endif

  BOOLPC
    (inputFromResults, m_input_from_results,
     "At the end of the simulation, write an input file containing the simulations "
     "parameters and generate a position file for restart purposes.");

/*m_properties.addProperty
    ("lambda", PropertyList::DOUBLE, &m_lambda_temp,
     new PLCDoubleRange(0, 1),
     "'lambda' value for modified velocity verlet algorithm.");*/

BOOLPC
    (randomize, m_randomize,
     "Initalize random number generator with a non-default seed. This makes the outcome of two simulations with identical setup different. All other modules that don't have their own randomize option, use randomize from Simulation.");
/*  BOOLPC
      (oneProp, m_one_prop,
       "Should the computed particle properties be common and unique for all pairs of "
           "species?\n"
           "Example: If we have two species A, B and 'oneProp' is set to 'true', we don't have "
           "separate local densities for (A, A), (B, B) and (A, B), but a SINGLE one, treating "
           "all particles identically.");*/
  m_properties.addProperty
    ("wtDistEps", PropertyList::DOUBLE, &c_wt_dist_eps, NULL,
     "(WallTriangle:) This is the distance a particle can have from a triangle "
     "in the plane of the triangle and still be detected as a hit. DO NOT CHANGE.");

  // the following is removed as attribute. The constant is set in
  // wall_triangle.cpp to '0.', which should be safe, if rsDispEps
  // (see below) is set to a value > 0

  /*  m_properties.addProperty
    ("wtTimeEps", PropertyList::DOUBLE, &c_wt_time_eps, NULL,
     "(WallTriangle:) The time to collision above which a collision is detected. "
     "DO NOT CHANGE.");*/


  m_properties.addProperty
    ("rsDispEps", PropertyList::DOUBLE, &c_rs_disp_eps, NULL,
     "(ReflectorStochastic:) The distance a particle is displaced from the wall "
     "when being reflected. DO NOT CHANGE.");
  m_properties.addProperty
    ("geomEps", PropertyList::DOUBLE, &g_geom_eps, NULL,
     "(Geometry) General epsilon for computational geometry. "
     "DO NOT CHANGE.");

    INTPC
    (stageIterations, m_stageSteps, 0,
     "Maximum number of iterations for runtime determination of Symbol stages.");

  m_name = "default";
  m_input_from_results = false;
//   m_lambda_temp = 0.5;
  m_randomize = false;
  m_stageSteps = 20;
//   m_one_prop = false;
#ifdef _OPENMP
global::n_threads = omp_get_num_procs();
#else
global::n_threads = 1;
#endif
}



void Simulation::setup()
{
  maxCutoff = -1;

  if (m_randomize) {
    m_rng.setSeed(getpid());
    MSG_DEBUG("Simulation::setup", "randomizing");
  }
  else {
    m_rng.setSeed(RNG_DEFAULT_SEED);
    MSG_DEBUG("Simulation::setup", "NOT randomizing");
  }

  NodeManyChildren::setup();

  if(!m_phase)
    throw gError("Simulation::setup", "No Phase defined.");

  m_particle_forces.resize(m_phase->manager()->nColours());
  for (size_t col = 0; col < m_phase->manager()->nColours(); ++col) {
    m_particle_forces[col] = new vector<GenF*>();
  }

  for (vector<GenF*>::iterator force = m_forces.begin(); force != m_forces.end(); force++) {
    if ((*force)->isParticleForce()) {
        m_particle_forces[(*force)->fCol()]->push_back(*force);
    }
    else if (!((*force)->isPairForce())) {
      m_other_forces.push_back(*force);
    }
  }

#ifdef _OPENMP
//   add copy-slot attributes
  for (size_t t = 0; t < global::n_threads; t++) {
    for (size_t col = 0; col < m_phase->manager()->nColours(); ++col) {
      Particle::s_tag_format[col].addAttribute
	("copy" + ObjToString(col) + ObjToString(t), 
	 DataFormat::VECTOR_DOUBLE,
	 /*persistent?*/true,
	 "ct");
    }
  }
#endif

}


void Simulation::writeWhenDone()
{
  string outName(m_name + "_out");
  string fileName(outName + ".in");
  string help;
  /* with the flags, it does not create output files (with gcc 3.x) */
  ofstream outFile(fileName.c_str());

  //    assert(!outFile.is_open());
  //    outFile.open(fileName.c_str()/*, ios::app | ios::binary | ios::trunc*/);
  if(!outFile)
    throw gError("Simulation::write: Unable to open file" + fileName);
  outFile.precision(8);

  NodeManyChildren::write(outFile);

  cout << "Input for new simulation written to file " << fileName << endl;
  outFile.close();
  //	return outFile;
}


Node *Simulation::instantiateChild(const string &name)
{
  Node *node;

  MSG_DEBUG("Simulation::instantiateChild", name);

//   MSG_DEBUG("Simulation::instantiateChild", "Valgrindcheck *m_controller: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(*m_controller));

//   MSG_DEBUG("Simulation::instantiateChild", "Valgrindcheck *m_phase: "<< VALGRIND_CHECK_VALUE_IS_DEFINED(*m_phase));

  if (name == "Phase") {
    if (m_phase)
      throw gError("Simulation::instantiateChild", "Only one Phase is allowed.");

/*        if (m_forces.empty())
          throw gError("Simulation::instantiateChild",
          "Phase must be defined after the forces");*/

    node = m_phase = new Phase(this);

    /* If more than one force, the max cutoff has to be given as argument
       and will be checked by the Phase... Ehm it's not exactly an argument
       to read, it's a public member of Simulation now, that's ugly so fixme!!! */
  } else if (name == "Controller") {
    if (m_controller)
      throw gError("Simulation::instantiateChild", "Only one Controller can be defined.");
    node = m_controller = new Controller(this);
  } else if (GenFType::exists(name)) {
    /* Is it a force? */
    GenF *force;

    //        MSG_DEBUG("Simulation::instantiateChild", "It's a force!");

    node = force = GenFType::byName(name).makeGenF(this);

    if (m_phase)
      throw gError
        ("Simulation::instantiateChild",
         "Forces must be defined before the Phase.");

    if (!m_controller)
      throw gError
        ("Simulation::instantiateChild",
         "Forces must be defined after the Controller.");

    // forces need the integrator as argument to be able to request a pointer
    // to its forces arrays


    m_forces.push_back(force);

  } else if (Meter_Factory::exists(name)) {
    /* Is it a meter? */
    Meter *meter;

    node = meter = Meter_Factory::byName(name).instantiate(this);

    if (!m_phase)
      throw gError
        ("Simulation::instantiateChild", name +
         " without Phase found. Meters have to be defined "
         "after the Phase they will measure.");
    m_meters.push_back(meter);
  } else if (Callable_Factory::exists(name)) {
    /* Is it a thermostat? */
    Callable *callable;

    if (!m_controller)
      throw gError
        ("Simulation::instantiateChild",
         "Callables must be defined after the Controller.");

    if (m_phase)
      throw gError
        ("Simulation::instantiateChild",
         "Callables must be defined before the Phase.");

    node = callable = Callable_Factory::byName(name).instantiate(this);

    m_callables.push_back(callable);
  } else if (WeightingFunction_Factory::exists(name)) {
    WeightingFunction *wf;

    node = wf = WeightingFunction_Factory::byName(name).instantiate(this);

    m_weighting_functions.push_back(wf);
  }
    else if (SymbolFactory::exists(name)) {
      if (!m_controller)
        throw gError
            ("Simulation::instantiateChild",
             "Symbols must be defined after the Controller.");

      if (m_phase)
	throw gError
		("Simulation::instantiateChild",
		"Symbols must be defined before the Phase.");

      Symbol *symbol;

      node = symbol = SymbolFactory::byName(name).instantiate(this);

      m_symbols.push_back(symbol);
  }
    else
    throw gError
      ("Simulation::instantiateChild", "'" + name + "' not found.");

  return node;
}


WeightingFunction *Simulation::findWeightingFunction(const string &name)
{
  FOR_EACH
    (vector<WeightingFunction*>,
     m_weighting_functions,
     if ((*__iFE)->name() == name)
       return *__iFE;
     );

  throw gError
    ("Simulation::findWeightingFunction",
     "No such weighting function defined: '" + name + "'");

  return NULL;
}

GenF *Simulation::findForceFunction(const string &name)
{
  vector<ColourPair*>::iterator __end = m_phase->manager()->colourPairs().end();
  for(vector<ColourPair*>::iterator __cp = m_phase->manager()->colourPairs().begin(); __cp != __end; ++__cp) {
    FOR_EACH
      (vector<GenF*>,
      (*(*__cp)->pairForces()),
      if ((*__iFE)->name() == name)
        return *__iFE;
      );
  }
  for (size_t col = 0; col < m_phase->manager()->nColours(); ++col) {
    FOR_EACH
      (vector<GenF*>,
      (*m_particle_forces[col]),
      if ((*__iFE)->name() == name)
        return *__iFE;
      );
  }

  FOR_EACH
    (vector<GenF*>,
     m_other_forces,
     if ((*__iFE)->name() == name)
       return *__iFE;
     );

  throw gError
    ("Simulation::findForceFunction",
     "No such force defined: '" + name + "'");

  return NULL;
}

GenF *Simulation::searchForceWithClassName(const string &name, const string& className)
{
  GenF* tmp = NULL;
  FOR_EACH
    (vector<GenF*>,
     m_forces,
     if ((*__iFE)->name() == name && (*__iFE)->className() == className) {
       if(tmp)   
	 throw gError
	   ("Simulation::findForceWithClassName:"+FILE_INFO,
	    "Two forces with same name \"" + name + "\" and class name \"" + className + "\" not allowed!");

       tmp = (*__iFE);
     }
     );

  return tmp;
}


