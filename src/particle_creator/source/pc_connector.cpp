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




#include <iostream>
#include <fstream>
#include "pc_connector.h"
// #include "gen_f.h"
// #include "gen_connector.h"
#include "f_specific.h"
// #include "gen_triplet.h"
#include "phase.h"
#include "colour_pair.h"
#include "simulation.h"
#include "manager_cell.h"
#include "triplet.h"
#include "quintet.h"

using namespace std;
#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation*) M_PHASE->parent())

/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleConnectorFile>
		particle_connector_file("ParticleConnectorFile");

ParticleConnectorFile::ParticleConnectorFile() :
	ParticleCreator() {
}

ParticleConnectorFile::ParticleConnectorFile(Boundary* boundary) :
	ParticleCreator(boundary) {
	init();
}

ParticleConnectorFile::~ParticleConnectorFile() {
}

//--- Methods ---
Particle* ParticleConnectorFile::getParticleFromNumber(size_t n,
		string freeOrFrozen) {
	vector<ParticleList*>::iterator it;
	vector<ParticleList*>::iterator end = m_particlesToConnect.end();
	size_t s, n0 = n;
	if (freeOrFrozen == "free") {
		for (it = m_particlesToConnect.begin(); it != end; ++it) {
			s = (**it).size();
			if (s > n) // contains n
				return &((**it)[n]);
			else {
				n -= s;
			}
		}
		throw gError("ParticleConnectorFile::getParticleFromNumber::"+FILE_INFO,
				string("Particle with index ") + ObjToString(n0) + " not found in lists. Did you define all the required species? Did you create the particles you refer to?");
	} else if (freeOrFrozen == "frozen") {
		end = m_frozenParticlesToConnect.end();
		for (it = m_frozenParticlesToConnect.begin(); it != end; ++it) {
			s = (**it).size();
			if (s > n) // contains n
				return &((**it)[n]);
			else {
				n -= s;
			}
		}
		throw gError("ParticleConnectorFile::getParticleFromNumber::"+FILE_INFO,
				string("Frozen particle with index ") + ObjToString(n0) + " not found in lists. Did you define all the required species? Did you create the particles you refer to?");
	} else
		throw gError("ParticleConnectorFile::getParticleFromNumber::", "Is the particle "
				+ ObjToString(n0) + " free or frozen?");
}

void ParticleConnectorFile::init() {
	m_properties.setClassName("ParticleConnectorFile");
	m_properties.setDescription("Reads information about particle connections from a file. The expected file format looks as in the following example:\n\n"
"pair spring solida solidb\n"
"triplet angularSpring\n"
"!!!\n"
"pair spring 120free 121free\n"
"triplet angularSpring 13free 14free 15free\n\n"
"The file is subdivided into two parts, the declaration of the required lists of particle-connections before the line \"!!!\" and the definition of the particle-connections themselves afterwards.\n"
"In the declaration part, pair and triplet lists must be declared.\nThe format for the pair lists is as follows: first, the identifier \"pair\", then the identifier of the list, then the two particle species. The format for the triplet lists is the same except that the particle species are NOT stated.\n"
"The format for the definition of the connections is as follows:\n"
"The first column indicates the type of connection. Currently allowed are \"pair\" and \"triplet\". The second column gives the identifier of the list, where the connection should be stored. Then the particle indices together with a statement whether it is a free or frozen particle follow. " 
"The particles must have been created by another ParticleCreator. "
		"NOTE: The particle indices start at 0 with the first "
		"species and continue with counting for the following ones! They do not restart at 0. BUT, the counting is performed SEPARATELY for free and frozen particles, both starting at 0!") ;

	//	STRINGPC(species, m_species, "Defines the species the ParticleCreator should connect."
	//			" It is allowed to specify multiple species separated by whitespace, e.g., "
	//			"\"H wall F\". The particle indices start at 0 with the first"
	//			"species and continue for the following ones.")
	//	;
	//	m_species = "UNDEF";

	STRINGPCINF
	(name, m_filename,
			"File containing the connector and specific force data. The file consists of one line"
			"per connector/force, where the first column is the force name and the "
			"following columns are the particle indices and 'frozen' or 'free' property.")
	;
	m_filename = "default.con";
}

void ParticleConnectorFile::setup() {

  s_createFrozenParts = false;
  
  // create the connected lists here;
  ManagerCell* manager = M_PHASE->manager();
  string type, species1, species2;
  // open file
  ifstream pos(m_filename.c_str());
  if (!pos)
    throw gError("ParticleConnectorFile::setup"+FILE_INFO,
		 "Error opening file '" + m_filename + "'.");  
  bool inDeclaration = true;
  while (inDeclaration) {
    pos >> skipws >> type;
    if(type == "pair") {
      pos >> skipws >> type;
      pos >> skipws >> species1;
      pos >> skipws >> species2;

//       MSG_DEBUG("ParticleConnectorFile::setup", "pair case: list=" << type << ", species1 = " << species1 << ", species2=" << species2);
      
      try {
ColourPair* cp = 
	M_PHASE->manager()->cp(manager->getColour(species1), manager->getColour(species2));
      /*size_t listIndex =*/ cp->createConnectedListIndex(type);
      }
      catch(gError& err) {
	throw gError("ParticleConnectorFile::setup", "Error while reading a line starting with \"pair\" in what is interpreted as the declaration part of " + m_filename + ". Does it have a declaration part at all? See the help text of ParticleConnectorFile. The error message was: " + err.message());
      }
    }
    else if(type == "triplet") {
      pos >> skipws >> type;
//       MSG_DEBUG("ParticleConnectorFile::setup", "found triplet, creating list " << type);
      try {
	M_PHASE -> createTripletListIndex(type);
      }
      catch(gError& err) {
	throw gError("ParticleConnectorFile::setup", "Error while reading a line starting with \"triplet\" in the declarations part of " + m_filename + ". Does it have a declaration at all? See the help text of ParticleConnectorFile. The error message was: " + err.message());
      }

    }
 else if(type == "quintet") {
      pos >> skipws >> type;
//       MSG_DEBUG("ParticleConnectorFile::setup", "found quintet, creating list " << type);
      try {
	M_PHASE -> createQuintetListIndex(type);
      }
      catch(gError& err) {
	throw gError("ParticleConnectorFile::setup", "Error while reading a line starting with \"quintet\" in the declarations part of " + m_filename + ". Does it have a declaration at all? See the help text of ParticleConnectorFile. The error message was: " + err.message());
      }

    }
    else if(type == "!!!") {
      inDeclaration = false;
    }
    else 
      throw gError("ParticleConnectorFile::setup", "unknown type" + type + "in declaration part of file " + m_filename + "! Allowed expressions: \"pair\",\"triplet\" ,\"quintet\", and \"!!!\".");

  }
  // the rest of the file is skipped and read in setupAfterParticleCreation
  pos.close();
    
  // the connections must be created after the creation of the particles
  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);
}

void ParticleConnectorFile::setupAfterParticleCreation() {

  size_t np1, np2, np3, np4, np5 ;
  Particle *p1, *p2, *p3, *p4, *pc ;

  size_t c;
  string s, type;
  GenF* f;
  Phase *phase= M_PHASE;
  ManagerCell* manager = phase->manager();
  string freeOrFrozen = "free";
  
  // which species should we connect?
  stringstream ss(m_species, ios_base::in);
  while (ss >> skipws >> s) {
    MSG_DEBUG("ParticleConnectorFile::setupAfterParticleCreation",
	      "Connecting particles of species "+s);
    c = manager->getColour(s);
    m_particlesToConnect.push_back(&(phase->particles(c)));
    m_frozenParticlesToConnect.push_back(&(phase->frozenParticles(c)));
  }
  
  // open file
  MSG_DEBUG("ParticleConnectorFile::setupAfterParticleCreation"+FILE_INFO, "trying to open file " << m_filename << " ...");
  ifstream pos(m_filename.c_str());
  if (!pos)
    throw gError("ParticleConnectorFile::setupAfterParticleCreation"+FILE_INFO,
		 "Error opening file '" + m_filename + "'.");

  // skip declaration of connected lists which was already read in setup();
  bool inDeclaration = true;
  while (inDeclaration) {
    pos >> skipws >> type;
    if(type == "pair") {
      // do nothing, skip listname, species1, species2
      pos >> skipws >> type;
      pos >> skipws >> type;
      pos >> skipws >> type;
    }
    else if(type == "triplet")
      // do nothing, skip listname
      pos >> skipws >> type;
    else if(type == "quintet")
      // do nothing, skip listname
      pos >> skipws >> type;
    else if(type == "!!!") 
      inDeclaration = false;
    else 
      throw gError("ParticleConnectorFile::setup", "unknown type" + type + "in declaration part of file " + m_filename + "! Allowed expressions: \"pair\",\"triplet\" ,\"quintet\", and \"!!!\".");

  }
  // read list type
  while ((pos >> skipws >> type) && !pos.eof()) {
    try {
      pos.exceptions(std::ios::failbit | std::ios::badbit | std::ios::eofbit);

      // read list name
      pos >> skipws >> s;
      // read first particle (of a total of one to three)
      pos >> skipws >> np1 >> skipws >> freeOrFrozen;
      p1 = getParticleFromNumber(np1,freeOrFrozen);

      if(type == "single") { // FSpecific
	// search s in forces (for FSpecific-case)
	f = M_SIMULATION->searchForceWithClassName(s, "FSpecific");
	static_cast<Fspecific*>(f)->addParticleToForce(p1);
//  	MSG_DEBUG("ParticleConnectorFile::setupAfterParticleCreation",
//  		  "Connected with FSpecific \"" << s << "\": " << p1->mySlot << " at " << p1->r);
      }
      else if(type == "pair") { // pair
	
	pos >> skipws >> np2 >> skipws >> freeOrFrozen;
	p2 = getParticleFromNumber(np2,freeOrFrozen);
	// find ColourPair, add pair, and (if necessary) new connection
	ColourPair* cp = phase->manager()->cp(p1->c, p2->c);
	size_t listIndex = cp->connectedListIndex(s);
	cp->addPairToConnection(p1, p2, listIndex);
//  	MSG_DEBUG("ParticleConnectorFile::setupAfterParticleCreation",
//  		  "Connected pair: " << p1->mySlot << " at " << p1->r <<
//  		  " and "<< p2->mySlot << " at " << p2->r << ", particle numbers from input = " + ObjToString(np1) + ", " + ObjToString(np2) + ", " + ObjToString(np3));
      }
      else if(type == "triplet") { // triplet
	pos >> skipws >> np2 >> skipws >> freeOrFrozen;
	p2 = getParticleFromNumber(np2,freeOrFrozen);
	pos >> skipws >> np3 >> skipws >> freeOrFrozen;
	p3 = getParticleFromNumber(np3,freeOrFrozen);
	size_t listIndex = phase -> tripletListIndex(s);

	// safety check for CONVENTION6:
	// check for consistency of the list's species
	tripletList* trl = M_PHASE->returnTripletList(s);
	tripletList::iterator firstTr = trl->begin();
	// already some stored triplet inside?	
	if(firstTr != trl->end()) {
	  // consistency-check
	  if(p1->c != firstTr->a->c ||
	     p2->c != firstTr->b->c ||
	     p3->c != firstTr->c->c 
	     )
	    throw gError("ParticleConnectorFile::setupAfterParticleCreation", "There seems to be an inconsistency in the order of the species of the connected particles in list \"" + s + "\". This is currently not allowed. p1 =  " + ObjToString(p1->mySlot) + ", c1 = "+ ObjToString(p1->c) + ", p2 = " + ObjToString(p2->mySlot) + ", c2 = "+ ObjToString(p2->c)+ ", p3 = " + ObjToString(p3->mySlot) + ", c3 = "+ ObjToString(p3->c) + ", expected colours = " + ObjToString(firstTr->a->c) + ObjToString(firstTr->b->c) + ObjToString(firstTr->c->c) + ", particle numbers from input = " + ObjToString(np1) + ", " + ObjToString(np2) + ", " + ObjToString(np3));
	}

	phase->addTriplet(p1, p2, p3, listIndex);
//  	MSG_DEBUG("ParticleConnectorFile::setupAfterParticleCreation",
//  		  "Connected triplet: " << p1->mySlot << " at " << p1->r<<
//  		  " and "<< p2->mySlot << " at " << p2->r<<
//  		  " and "<< p3->mySlot << " at " << p3->r << ", colours: " << p1->c << p2->c << p3->c << ", slots = (" << p1->mySlot << ", " << p2->mySlot << ", " << p3->mySlot << ")" << ", particle numbers from input = " + ObjToString(np1) + ", " + ObjToString(np2) + ", " + ObjToString(np3));
	
      }
      else if(type == "quintet") { // quintet
	pos >> skipws >> np2 >> skipws >> freeOrFrozen;
	p2 = getParticleFromNumber(np2,freeOrFrozen);
	pos >> skipws >> np3 >> skipws >> freeOrFrozen;
	p3 = getParticleFromNumber(np3,freeOrFrozen);
	pos >> skipws >> np4 >> skipws >> freeOrFrozen;
	p4 = getParticleFromNumber(np4,freeOrFrozen);
	pos >> skipws >> np5 >> skipws >> freeOrFrozen;
	pc = getParticleFromNumber(np5,freeOrFrozen);
	size_t listIndex = phase -> quintetListIndex(s);

	// safety check for CONVENTION6:
	// check for consistency of the list's species
	quintetList* quin = M_PHASE->returnQuintetList(s);
	quintetList::iterator firstQuin = quin->begin();
	// already some stored quintet inside?	
	if(firstQuin != quin->end()) {
	  // consistency-check
	  if(p1->c != firstQuin->p00->c ||
	     p2->c != firstQuin->p20->c ||
    	     p3->c != firstQuin->p22->c ||
	     p4->c != firstQuin->p02->c ||
	     pc->c != firstQuin->p11->c 
	     )
	    throw gError("ParticleConnectorFile::setupAfterParticleCreation", "There seems to be an inconsistency in the order of the species of the connected particles in list \"" + s + "\". This is currently not allowed. p1 =  " + ObjToString(p1->mySlot) + ", c1 = "+ ObjToString(p1->c) + ", p2 = " + ObjToString(p2->mySlot) + ", c2 = "+ ObjToString(p2->c)+ ", p3 = " + ObjToString(p3->mySlot) + ", c3 = "+ ObjToString(p3->c) + ", p4 = " + ObjToString(p4->mySlot) + ", c4 = "+ ObjToString(p4->c)+ ", pc = " + ObjToString(pc->mySlot) + ", cc = "+ ObjToString(pc->c)+ ", expected colours = " + ObjToString(firstQuin->p00->c) + ObjToString(firstQuin->p20->c) + ObjToString(firstQuin->p22->c) + ObjToString(firstQuin->p02->c) + ObjToString(firstQuin->p11->c)+ ", particle numbers from input = " + ObjToString(np1) + ", " + ObjToString(np2) + ", " + ObjToString(np3)+ ", " + ObjToString(np4)+ ", " + ObjToString(np5));
	}

	phase->addQuintet(p1, p2, p3, p4, pc, listIndex);	
      }
      else // nothing found
	throw gError("ParticleConnectorFile::setupAfterParticleCreation"+FILE_INFO,
		     "Cannot handle connection type \"" + type + "\"!");
      pos.exceptions(std::ios::goodbit);
      
    }
    catch(ifstream::failure e) {
      throw gError("ParticleConnectorFile::setupAfterParticleCreation::"+FILE_INFO,
		   "Error reading file '" + m_filename + "': " + e.what());
      
    }
  } // end of while
  
  pos.close();

}
