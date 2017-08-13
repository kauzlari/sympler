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

#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include "cell.h"
#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "particle_cache.h"
#include "boundary.h"
#include "pc_file.h"
using namespace std;
/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorFile>
		particle_creator_file("ParticleCreatorFile");

#define M_BOUNDARY ((Boundary *) m_parent)
#define M_PHASE ((Phase *) M_BOUNDARY->parent())
#define M_MANAGER M_PHASE->manager()

//---- Constructor/Destructor ----

ParticleCreatorFile::ParticleCreatorFile(Boundary *boundary) :
	ParticleCreatorFreePCalc(boundary) {
	init();
}

//---- Methods ----

string ParticleCreatorFile::readNext(ifstream &pos) {
  int level = 0;
  char c;
  
  string s("");
  
  pos >> skipws;
  
  while ((c = pos.get()) == ' ');
  
  while (level > 0 || isalnum(c) || c == '-' || c == '.' || c == '(' || c == ')') {
//     MSG_DEBUG("ParticleCreatorFile::readNext(ifstream &pos)", "c = " << c);
    if (c == '(')
      level++;
    else if (c == ')')
      level--;
    s += c;
    
    c = pos.get();
  }
  
  /*
    if (pos.peek() == '(') {
    char c;
    while ((c = pos.get()) != ')')
    s += c;
    s += ")";
    } else
    pos >> s;
  */
//    MSG_DEBUG("ParticleCreatorFile::readNext(ifstream &pos)", "END: c = " << c << "; s = " << s);
  
  return s;
}

void ParticleCreatorFile::createParticles(){
  MSG_DEBUG("ParticleCreatorFile::createParticles","start");
    
  Phase *phase = M_PHASE;
  ManagerCell *manager = M_MANAGER;
  ifstream pos(m_filename.c_str());
  string s;
  vector<list<pair<string, int> > > tags;
  vector<list<bool> > writeTags;
  string species;

  initTransform();

  if (!pos)
    throw gError("ParticleCreatorFile::createParticles", "Error opening file '" + m_filename + "'.");

  tags.resize(manager->nColours());
  writeTags.resize(manager->nColours());

  pos >> skipws >> s;
  // 1st while
  while (s != "!!!" && !pos.eof()) {
MSG_DEBUG("ParticleCreatorFile::createParticles","start of 1st file-read-while; s = " << s);

    size_t c;

    species = s;
    if (m_species == "UNDEF" || m_species == species){ 
      c = manager->getColour(s);
      MSG_DEBUG("ParticleCreatorFile::createParticles", "requesting colour for species" << s << " and got " << c);
    }

    pos >> skipws >> s;
//     MSG_DEBUG("ParticleCreatorFile::createParticles","before 2nd file-read-while; s = " << s);

    // 2nd while: this loop is just entered if there are declared tag-fields
    while (s != "!!!" && !pos.eof()) 
    {
MSG_DEBUG("ParticleCreatorFile::createParticles","start of 2nd file-read-while");
      if (m_species == "UNDEF" || m_species == species) { 
MSG_DEBUG("ParticleCreatorFile::createParticles","in if " << species);
	
	size_t slot, offset;
        bool tempBool = true;

	MSG_DEBUG("ParticleCreatorFile::createParticles", "exists-check: now at string " << s);
	if(Particle::s_tag_format[c].attrExists(s)) {
	  slot = Particle::s_tag_format[c].attrByName(s).index;
	  offset = Particle::s_tag_format[c].attrByName(s).offset;
	  MSG_DEBUG("ParticleCreatorFile::createParticles", "EXISTS: now at string " << s << " with offset " << offset << " and slot(index) " << slot);
	}
	// don't write because not used in this simulation
	else {
	  tempBool = false;
	  MSG_INFO("ParticleCreatorFile::createParticles", "Skipping column \"" << s << "\" from particle-file because apparently not used in the current simulation.");
	}
	
	// add in any case, so that we can later on decide whether to skip or to write
	// c always has a meaningful value because of if (m_species == "UNDEF" || m_species == species)
        tags[c].push_back(pair<string, int>(s, slot));
        
//         MSG_DEBUG
//           ("ParticleCreatorFile::createParticles",
//           "Found tag '" << s << "' for species '" << species << "'. Now"
//           " searching for ValCalculator for this tag");
        
	if(tempBool) {
	  
	  pair<size_t, size_t> theSlots;
	  // in the following loops we check, whether there exists a
	  // ValCalculator for the found tag. If yes, this ParticleCreator 
	  // WILL NOT write the values from the file into the particles, since 
	  // the ValCalculator will compute them during update of the PairList

	  // first, loop over ColourPairs
	  FOR_EACH_COLOUR_PAIR
	    (
	     manager,

	     // next "if" not needed anymore (2011-05-18) because Calcs from any colour may generally write in any c 
// 	     if(cp->firstColour() == c || cp->secondColour() == c)
// 	       {
		 // MSG_DEBUG("ParticleCreatorFile::createParticles", "CP correct" << cp->toString());
		 
		 // loop over stages: should be OK, because Symbols are already sorted
		 // see in Controller::run()
		 
		 for(size_t stage = 0; (stage <= cp->maxStage()/*ColourPair::s_maxStage*/ /*VC_MAX_STAGE*/ && tempBool); ++stage)
		   {
		     // MSG_DEBUG("ParticleCreatorFile::createParticles", "stage = " << stage);
		     
		     vector<ValCalculator*>& vCs = cp->valCalculators(stage);
#ifdef _OPENMP
		     vector<ValCalculator*> vCPs = cp->valCalculatorParts(stage);
#endif
		     // loop over ValCalculators
		     for(vector<ValCalculator*>::iterator vCIt = vCs.begin(); 
			 (vCIt != vCs.end() && tempBool); ++vCIt)
		       {
			 // MSG_DEBUG("ParticleCreatorFile::createParticles", "now VC = " << (*vCIt)->myName());
			 // the FOR_EACH_COLOUR_PAIR does not seem to like the comma 
			 // in pair<...
			 //               pair<size_t, size_t> theSlots;
			 try
			   {
			     (*vCIt)->mySlots(theSlots);
			     // MSG_DEBUG("ParticleCreatorFile::createParticles", "theSlots = (" << theSlots.first << ", " << theSlots.second << "), offset = " << offset );


			     // next "if" not needed anymore (2011-05-18) because Calcs from any cp may generally write in any c 
// 			     if(cp->firstColour() == c)
			     if(manager->getColour((*vCIt)->firstWriteSpecies()) == c)
			       {
				 if(theSlots.first == offset) 
				   {
				     if((*vCIt)->mySymbolName() != s) {
				       MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol \"" + s + "\" from particle-file! Found ValCalculator calculating symbol \"" + (*vCIt)->mySymbolName() + "\" at same memory position. Contact the programmers. Aborting.");
				       abort();
				     }
				     if(!(*vCIt)->doesOverwrite()) {
				       tempBool = false;
				       MSG_INFO("ParticleCreatorFile::createParticles", 
						 "found " << (*vCIt)->myName() 
						 << " for symbol " << (*vCIt)->mySymbolName() 
						 << ". Skipping the corresponding column-values in" 
						 " the particle file.");
				     }
				   }
			       }
			     else
			       {
			     // next "if" now needed (2011-05-18) because Calcs from any cp may generally write in any c 
				 if (manager->getColour((*vCIt)->secondWriteSpecies()) == c) {

				   if(theSlots.second == offset) {
				     if((*vCIt)->mySymbolName() != s) {
				       MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol \"" + s + "\" from particle-file! Found ValCalculator calculating symbol \"" + (*vCIt)->mySymbolName() + "\" at same memory position. Contact the programmers. Aborting.");
				       abort();
				     }
				     if(!(*vCIt)->doesOverwrite()) {
				       tempBool = false;
				       MSG_INFO("ParticleCreatorFile::createParticles", 
						 "found " << (*vCIt)->myName() 
						 << " for symbol " << (*vCIt)->mySymbolName() 
						 << ". Skipping the corresponding column-values in" 
						 " the particle file.");               
				     }
				   }
				 } // end of if (manager->getColour((*vCIt)->secondWriteSpecies()) == c)
			       } // end of else of if ((*vCIt)->manager->getColour(firstWriteSpecies()) == c)
			   } // end of try
			 // that's if (*vCIt)->mySlots(theSlots) throws an exception
			 catch(gError& err)
			   {
// 			     MSG_DEBUG("ParticleCreatorFile::createParticles", 
// 				       "Ignoring ValCalculatorPair: " 
// 				       << (*vCIt)->mySymbolName());
			   }                
		       } // end of for(vector<ValCalculator*>::iterator vCIt...
#ifdef _OPENMP
		     for(vector<ValCalculator*>::iterator vCIt = vCPs.begin(); 
			 (vCIt != vCPs.end() && tempBool); ++vCIt) {
		       // MSG_DEBUG("ParticleCreatorFile::createParticles", "now VC = " << (*vCIt)->myName());
		       // the FOR_EACH_COLOUR_PAIR does not seem to like the comma 
		       // in pair<...
		       //               pair<size_t, size_t> theSlots;
		       try {
			 (*vCIt)->mySlots(theSlots);
			 // MSG_DEBUG("ParticleCreatorFile::createParticles", "theSlots = (" << theSlots.first << ", " << theSlots.second << "), offset = " << offset );


			 // next "if" not needed anymore (2011-05-18) because Calcs from any cp may generally write in any c 
			 // if(cp->firstColour() == c) {
			 if(manager->getColour((*vCIt)->firstWriteSpecies()) == c) {

			   if(theSlots.first == offset) {
			     if((*vCIt)->mySymbolName() != s) {
			       MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ValCalculator calculating symbol " + (*vCIt)->mySymbolName() + " at same memory position. Contact the programmers.");
			       abort();
			     }
			     if(!(*vCIt)->doesOverwrite()) {
			       tempBool = false;
			       MSG_INFO("ParticleCreatorFile::createParticles", 
					 "found " << (*vCIt)->myName() 
					 << " for symbol " << (*vCIt)->mySymbolName() 
					 << ". Skipping the corresponding column-values in" 
					 " the particle file.");
			     }
			     
			   } // end of if(theSlots.first == offset)
			 } // end of if(manager->getColour((*vCIt)->firstWriteSpecies()) == c)
			 else {

			   // next "if" now needed (2011-05-18) because Calcs from any cp may generally write in any c 
			   if (manager->getColour((*vCIt)->secondWriteSpecies()) == c) {

			     if(theSlots.second == offset) {
			       if((*vCIt)->mySymbolName() != s) {
				 MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ValCalculator calculating symbol " + (*vCIt)->mySymbolName() + " at same memory position. Contact the programmers.");
				 abort();
			       }
			       if(!(*vCIt)->doesOverwrite()) {
				 tempBool = false;
				 MSG_INFO("ParticleCreatorFile::createParticles", 
					   "found " << (*vCIt)->myName() 
					   << " for symbol " << (*vCIt)->mySymbolName() 
					   << ". Skipping the corresponding column-values in" 
					   " the particle file.");               
			       }
			     }
			   }
			 } // end of else of if (manager->getColour((*vCIt)->firstWriteSpecies()) == c)
			 
		       } // end of try ...
		       catch(gError& err) {
// 			 MSG_DEBUG("ParticleCreatorFile::createParticles", 
// 				   "Ignoring ValCalculatorPair: " 
// 				   << (*vCIt)->mySymbolName());
		       }
		     } // end of for(vector<ValCalculator*>::iterator vCIt ...) ...
#endif
			 
		   } // end of for(size_t stage = 0 ...)... 
		 for(size_t stage = 0; (stage <= cp->maxStage_0()/*ColourPair::s_maxStage*/ /*VC_MAX_STAGE*/ && tempBool); ++stage) {
		   // MSG_DEBUG("ParticleCreatorFile::createParticles", "stage = " << stage);
		   
		   vector<ValCalculator*>& vCs = cp->valCalculators_0(stage);
#ifdef _OPENMP
		   vector<ValCalculator*>& vCPs = cp->valCalculatorParts_0(stage);
#endif
		   // loop over ValCalculators
		   for(vector<ValCalculator*>::iterator vCIt = vCs.begin(); 
		       (vCIt != vCs.end() && tempBool); ++vCIt) {
		     // MSG_DEBUG("ParticleCreatorFile::createParticles", "now VC = " << (*vCIt)->myName());
		     // the FOR_EACH_COLOUR_PAIR does not seem to like the comma 
		     // in pair<...
		     //               pair<size_t, size_t> theSlots;
		     try {
		       (*vCIt)->mySlots(theSlots);
		       // MSG_DEBUG("ParticleCreatorFile::createParticles", "theSlots = (" << theSlots.first << ", " << theSlots.second << "), offset = " << offset );




		       // next "if" not needed anymore (2011-05-18) because Calcs from any cp may generally write in any c 
		       // if(cp->firstColour() == c) {
		       if(manager->getColour((*vCIt)->firstWriteSpecies()) == c) {

			 if(theSlots.first == offset) {
			   if((*vCIt)->mySymbolName() != s) {
			     MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ValCalculator calculating symbol " + (*vCIt)->mySymbolName() + " at same memory position. Contact the programmers.");
			     abort();
			   }
			   if(!(*vCIt)->doesOverwrite()) {
			     tempBool = false;
			     MSG_INFO("ParticleCreatorFile::createParticles", 
				       "found " << (*vCIt)->myName() 
				       << " for symbol " << (*vCIt)->mySymbolName() 
				       << ". Skipping the corresponding column-values in" 
				       " the particle file.");
			   }
			 }
		       }
		       else {

			 // next "if" now needed (2011-05-18) because Calcs from any cp may generally write in any c 
			 if (manager->getColour((*vCIt)->secondWriteSpecies()) == c) {

			   if(theSlots.second == offset) {
			     if((*vCIt)->mySymbolName() != s) {
			       MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ValCalculator calculating symbol " + (*vCIt)->mySymbolName() + " at same memory position. Contact the programmers.");
			       abort();
			     }
			     if(!(*vCIt)->doesOverwrite()) {
			       tempBool = false;
			       MSG_INFO("ParticleCreatorFile::createParticles", 
					 "found " << (*vCIt)->myName() 
					 << " for symbol " << (*vCIt)->mySymbolName() 
					 << ". Skipping the corresponding column-values in" 
					 " the particle file.");               
			     }
			   }
			 }
		       }
		     }
		     catch(gError& err) {
// 		       MSG_DEBUG("ParticleCreatorFile::createParticles", 
// 				 "Ignoring ValCalculatorPair: " 
// 				 << (*vCIt)->mySymbolName());
		     }
		   }
#ifdef _OPENMP
		   for(vector<ValCalculator*>::iterator vCIt = vCPs.begin(); 
		       (vCIt != vCPs.end() && tempBool); ++vCIt) {
		     // MSG_DEBUG("ParticleCreatorFile::createParticles", "now VC = " << (*vCIt)->myName());
		     // the FOR_EACH_COLOUR_PAIR does not seem to like the comma 
		     // in pair<...
		     //               pair<size_t, size_t> theSlots;
		     try {
		       (*vCIt)->mySlots(theSlots);
		       // MSG_DEBUG("ParticleCreatorFile::createParticles", "theSlots = (" << theSlots.first << ", " << theSlots.second << "), offset = " << offset );

		       // next "if" not needed anymore (2011-05-18) because Calcs from any cp may generally write in any c 
		       // if(cp->firstColour() == c) {
		       if(manager->getColour((*vCIt)->firstWriteSpecies()) == c) {

			 if(theSlots.first == offset) {
			   if((*vCIt)->mySymbolName() != s) {
			     MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ValCalculator calculating symbol " + (*vCIt)->mySymbolName() + " at same memory position. Contact the programmers.");
			     abort();
			   }
			   if(!(*vCIt)->doesOverwrite()) {
			     tempBool = false;
			     MSG_INFO("ParticleCreatorFile::createParticles", 
				       "found " << (*vCIt)->myName() 
				       << " for symbol " << (*vCIt)->mySymbolName() 
				       << ". Skipping the corresponding column-values in" 
				       " the particle file.");
			   }
			 }
		       }
		       else {

			 // next "if" now needed (2011-05-18) because Calcs from any cp may generally write in any c 
			 if (manager->getColour((*vCIt)->secondWriteSpecies()) == c) {

			   if(theSlots.second == offset) {
			     if((*vCIt)->mySymbolName() != s) {
			       MSG_DEBUG("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ValCalculator calculating symbol " + (*vCIt)->mySymbolName() + " at same memory position. Contact the programmers.");
			       abort();
			     }
			     if(!(*vCIt)->doesOverwrite()) {
			       tempBool = false;
			       MSG_INFO("ParticleCreatorFile::createParticles", 
					 "found " << (*vCIt)->myName() 
					 << " for symbol " << (*vCIt)->mySymbolName() 
					 << ". Skipping the corresponding column-values in" 
					 " the particle file.");               
			     }
			   }
			 }
		       }
		     }
		     catch(gError& err) {
// 		       MSG_DEBUG("ParticleCreatorFile::createParticles", 
// 				 "Ignoring ValCalculatorPair: " 
// 				 << (*vCIt)->mySymbolName());
		     }
		   }
#endif
		 }
		     
		 // can we stop FOR_EACH_COLOUR_PAIR now?
		 if(!tempBool)
		   {
		     __cp = __end;
		     // important because there still comes the ++__cp from the loop
		     --__cp;
		   }
// 	       } // end if(cp->firstColour == c || cp->secondColour == c) (see above why commented out)
	     );
	     
	     // do we also have to search in the ParticleCaches? ...
	     if(tempBool)
	       {// ... yes we have to search in the ParticleCaches

		 // also loop over colours now (2011-05-18) because a ParticleCache my generally write in any colour;
		 for(size_t col = 0; col < phase->nColours(); ++col) {

		   
		   // loop over stages
		   for(size_t stage = 0; stage <= Particle::s_maxStage && tempBool; ++stage)
		     {
		       if (Particle::s_cached_properties[col].size() > stage) {
			 FOR_EACH
			   (vector<ParticleCache*>,
			    Particle::s_cached_properties[col][stage],

			    // also check colours now (2011-05-18) because a ParticleCache my generally write in any colour;
			    if ((*(*__iFE)).writeColour() == c && (*(*__iFE)).offset() == offset) {
				if((*(*__iFE)).mySymbolName() != s)
				  throw gError("ParticleCreatorFile::createParticles", "Fatal error for symbol " + s + "from particle-file! Found ParticleCache (stage 1) calculating symbol " + (*(*__iFE)).mySymbolName() + " at same memory position. Contact the programmers. Read-species of this ParticleCache: " + manager->species(col) + ". Write-species of this ParticleCache: " + manager->species(c));
				else {
				  if(!(*(*__iFE)).doesOverwrite()) {
				    MSG_INFO("ParticleCreatorFile::createParticles", 
					      "found " << (*(*__iFE)).myName() 
					      << " for symbol " << (*(*__iFE)).mySymbolName() 
					      << ". Skipping the corresponding column-values in" 
					      " the particle file.");
				    tempBool = false;
				  }
				}
				
			    }
			    // may we abort the loop over the ParticleCaches?
			    if(!tempBool)
			      {
				__iFE = __end;
				// important because there still comes the ++__iFE from the loop
				--__iFE; 
			      }
			    );
		       }
		     if (Particle::s_cached_properties_0[col].size() > stage) 
		       {
			 FOR_EACH
			   (
			    
			    vector<ParticleCache*>,
			    Particle::s_cached_properties_0[col][stage],

			    // also check colours now (2011-05-18) because a ParticleCache my generally write in any colour;
			    if ((*(*__iFE)).writeColour() == c && (*(*__iFE)).offset() == offset) {
			      if((*(*__iFE)).mySymbolName() != s)
				throw gError("ParticleCreatorFile::createParticles", "Fatal error for symbol \"" + s + "\" from particle-file! Found ParticleCache (stage 0) calculating symbol \"" + (*(*__iFE)).mySymbolName() + "\" at same memory position. Contact the programmers. Read-species of this ParticleCache: " + manager->species(col) + ". Write-species of this ParticleCache: " + manager->species(c));
			      else {
				  if(!(*(*__iFE)).doesOverwrite()) {
				    MSG_INFO("ParticleCreatorFile::createParticles", 
					      "found " << (*(*__iFE)).myName() 
					      << " for symbol " << (*(*__iFE)).mySymbolName() 
					      << ". Skipping the corresponding column-values in" 
					      " the particle file.");
				    tempBool = false;
				  }
			      }
			      
			    }
			    // may we abort the loop over the ParticleCaches?
			    if(!tempBool)
			      {
				__iFE = __end;
				// important because there still comes the ++__iFE from the loop
				--__iFE; 
			      }
			    );
		       }
		     } // end of for(size_t stage = 0...)
		 } // end of for(size_t col = 0; col < phase->nColours(); ++col)
		 
	       }
	} // end of if(tempBool) for search in Calculators and Caches if true
	if(tempBool) 
	  MSG_INFO("ParticleCreatorFile::createParticles", "I will write the tag " << s);

	// c always meaningful because of if (m_species == "UNDEF" || m_species == species)
	writeTags[c].push_back(tempBool);
	
      } // end of if (m_species == "UNDEF" || m_species == species)
      
      pos >> skipws >> s;

// MSG_DEBUG("ParticleCreatorFile::createParticles", "end of 2nd while (s != \"!!!\" && !pos.eof()); s = " << s);

    } // end of 2nd while (s != "!!!" && !pos.eof()) (just entered if tag-fields have been declared in the file)
    
    if (pos.eof())
      throw gError("ParticleCreatorFile::createParticles", "File corrupted.");
    
    pos >> skipws >> s;
    
  } // end of 1st while (s != "!!!" && !pos.eof())
  
  // initial read for species of real particle data starting now   
  pos >> skipws >> species;
  // 3rd while (for the real particle data)

  while (species != "!!!" && !pos.eof()) {

//     MSG_DEBUG("ParticleCreatorFile::createParticles", "start of 3rd while (species != \"!!!\" && !pos.eof()); species =  " << species << "; m_species = " << m_species);

    Particle p;
    Cell *c;
    size_t colour;
    string freeOrFrozen = "free";
    
    // if this is a species to be ignored (false case), read to end of line without doing anything
    //     bool reallyCreate;
    if (m_species == "UNDEF" || m_species == species) {
      colour = manager->getColour(species);
    
      
      readParticle(p,pos,freeOrFrozen);
  
      p.setColour(colour);
    

      // old (BUGGY! 2013-07-29) style of species checking
      //     if (m_species == "UNDEF" || m_species == species){
      list<bool>::iterator boolIt = writeTags[colour].begin();
      for (list<pair<string, int> >::iterator j = tags[colour].begin(); j != tags[colour].end(); j++) {

//  	MSG_DEBUG("ParticleCreatorFile::createParticles", "tags loop; colour = " << colour << "; tagname = " << j->first << "; s before readNext = " << s);

	s = readNext(pos);

// 	MSG_DEBUG("ParticleCreatorFile::createParticles", "before p.tag.fromStringByIndex: now at symbol " << j->first << ", index = " << j->second << ", just read s = " << s);
	
	if(*boolIt) {
	  p.tag.fromStringByIndex(j->second, s);
	}
	++boolIt;
      }
      
      transformPos(p);
      
      c = manager->findCell(p.r);
      if (c) {
    
        if(m_particlesInside){
         if (M_BOUNDARY->isInside(p.r)) {
           p.g = c->group();
	  
           if (freeOrFrozen == "frozen")
	    m_particles_frozen[p.g].newEntry() = p;
	  else
	    m_particles[p.g].newEntry() = p;
          }
	}
        if(!m_particlesInside){        
         if (!(M_BOUNDARY->isInside(p.r))) {
           p.g = c->group();
	  
           if (freeOrFrozen == "frozen")
	    m_particles_frozen[p.g].newEntry() = p;
	  else
	    m_particles[p.g].newEntry() = p;
          }
	}
      }
      //     } // end of if(m_species ...)
    } // end of if(m_species ...)
    else {// ignore this line since species not used in this simulation
      getline(pos, species);
//       MSG_DEBUG("ParticleCreatorFile::createParticles", "ignored the species and read the whole line '" << species << "' without doing anything.");
    }

    // read species at beginning of next line for next round of loop
    pos >> skipws >> species;
//     MSG_DEBUG("ParticleCreatorFile::createParticles", "end of 3rd while(species != \"!!!\" && !pos.eof()); next species = " << species << "; s = " << s << ", p.v.z = " << p.v.z);
  } // end of 3rd while(species != "!!!" ...) (for the real particle data)
  pos.close();
  /* next will call the one in ParticleCreatorFree */
  flushParticles();
}


void ParticleCreatorFile::readParticle(Particle &p, ifstream &pos, string &freeOrFrozen ){
    
      pos >> skipws >> freeOrFrozen;
      
      if (freeOrFrozen == "free" || freeOrFrozen == "frozen") {
	pos >> skipws >> p.r.x >> skipws >> p.r.y >> skipws >> p.r.z
	    >> skipws >> p.v.x >> skipws >> p.v.y >> skipws >> p.v.z;
      } else {
	p.r.x = atof(freeOrFrozen.c_str());
	pos >> skipws >> p.r.y >> skipws >> p.r.z >> skipws >> p.v.x
	    >> skipws >> p.v.y >> skipws >> p.v.z;
	freeOrFrozen = "free";
      }
        pos.ignore(HUGE_VAL,'\n');
    
 }

void ParticleCreatorFile::setup() {
  //	ParticleCreatorFree::setup();
  // this PC uses m_species as filter for the given input file
 
    myCutoff = M_PHASE->pairCreator()->interactionCutoff(); //as in pc_wall for adjustBoxSize()
      if (m_species != "UNDEF") {
    m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
  } else {
    //   m_colour = ALL_COLOURS;
    ifstream pos(m_filename.c_str());
    string s;
    
    size_t c;
    if (!pos)
      throw gError("ParticleCreatorFile::createParticles", "Error opening file '" + m_filename + "'.");
    
    pos >> skipws >> s;
    while (s != "!!!" && !pos.eof()) {
      
      c = M_MANAGER->getColour/*AndAdd*/(s);
      pos >> skipws >> s;
      
      while (s != "!!!" && !pos.eof())
	pos >> skipws >> s;
      pos >> skipws >> s;
      
    }
    pos.close();
    
}
}


void ParticleCreatorFile::init() {
  m_properties.setClassName("ParticleCreatorFile");
  
  m_properties.setDescription
    ("Loads a particle configuration from a file. The file format is the one obtained"
     " by setting "
     "'inputFromResults = yes' in the Simulation object.\n"
     "If the attribute 'species' is defined, ParticleCreatorFile uses it as a filter"
     " in order to extract all files of that 'species' from the given file."
     " Otherwise the particles from all species that have been defined in the XML-input"
     " will be extracted from the configuration file.\n"
     "Frozen or free particles can be created by stating 'frozen' or 'free' after the name of the species."
     " If nothing is stated, free particles will be created (in order to be backward compatible).\n"
     "Additional user defined Symbols can be given as well. Differently to positions or velocities"
     " they must be given as one comma-separated string inside of \"()\" brackets as shown in the"
     " example below for a vector. A 3x3 tensor is defined in the same way with additional brackets as"
     " a list of three vectors.\n "
     "Also note: If your file contains additional symbols to be assigned to the particles, it might be a good idea to check the initialisation output during the run for \"Skipping...\" status reports of ParticleCreatorFile. The module will skip assigning values of those columns of your file for which you have defined 'Symbol' modules with attribute 'overwrite = \"no\"' in the XML-input file.\n"
     "The module will skip those columns COMPLETELY for which the symbol does not seem to be used at all in the XML-input. This means, declarations of additional degrees of freedom still have to be performed in the XML-input.\n"
     "Example file:\n"
     "--- BEGINNING OF FILE ---\n"
     "H2O q0 q1 q2 q3 omega !!!\n"
     "D2O q0 q1 q2 q3 omega !!!\n"
     "!!!\n"
     "H2O free 1.385 1.385 4.155 2.7705274 1.7009424 2.8508353 0.707406 -0.139236 -0.0898079 0.687113 (-0.0511271, -0.107274, -0.0166022)\n"
     "H2O free 1.385 1.385 6.925 1.4525585 1.1747588 -1.2542176 0.74004 -0.122624 -0.0572537 0.658807 (-0.0185539, -0.0782408, 0.0170039)\n"
     "D2O free 26.315 26.315 23.545 1.1201054 0.37740195 -1.1816764 0.997009 0.0433965 0.0465469 0.0438526 (-0.0561241, 0.167672, 0.012043)\n"
     "D2O free 26.315 26.315 26.315 -2.2959708 2.5165383 2.3576235 0.953161 0.262601 -0.0545917 -0.139802 (-0.0188484, -0.0760585, 0.070389)\n"
     "!!!\n"
     "--- END OF FILE ---\n"
     "First, each new line introduces one particle species by giving the species name, and then followed by the names of additional attributes, all space separated and terminated by '!!!'.\n"
     "Then, in another new line this section is terminated by another '!!!'.\n"
     "The particles are defined one per row, starting with their species, optionally followed by the label \"free\" or \"frozen\", and then followed by three position-values, three velocity values, and then the values of the additional attributes in the order specified in the header.\n"
     "After the last particle, the file is terminated by another new line containing '!!!'.\n"
  
     );

	STRINGPCINF
	(name, m_filename,
			"File containing the position and velocity information.")
	;
        BOOLPC
                (particlesInside, m_particlesInside,
     " true if particles are inside inside of geometry (fluid particles), false"
                " if they are outside (wall particles). Useful when working with STL geometries. Default is true so that it is backwards compatible. ");
        
	m_filename = "default.pos";
        m_particlesInside = "true";
      
}

void ParticleCreatorFile::flushParticles() {
	Phase *phase= M_PHASE;
	size_t counter = 0;
	if (!m_particles_frozen.empty()) {
	  for (map<int, ParticleList>::iterator g = m_particles_frozen.begin(); g
		 != m_particles_frozen.end(); g++) {
	    SL_FOR_EACH
	      (Particle, g->second,
	       transformVel(*__iSLFE);
	       phase->addFrozenParticle(*__iSLFE);
	       ++counter;
	       );
	  }
	}

	if (!m_particles.empty()) {
	  for (map<int, ParticleList>::iterator g = m_particles.begin(); g
		 != m_particles.end(); g++) {
	    SL_FOR_EACH
	      (Particle, g->second,
	       transformVel(*__iSLFE);
	       phase->addParticle(*__iSLFE);
	       ++counter;
	       )
	      ;
	  }
	}
	MSG_DEBUG("ParticleCreatorFile::flushParticles" , counter << " particles added");
	m_particles_frozen.clear();
	m_particles.clear();
}
void ParticleCreatorFile::adjustBoxSize(point_t &size, bool_point_t& frameRCfront,  bool_point_t& frameRCend){

    ifstream pos(m_filename.c_str());
    string s ;
    pos >> skipws >> s;
        while (s != "!!!" && !pos.eof()) {

          pos >> skipws >> s;
          while (s != "!!!" && !pos.eof()){
            pos >> skipws >> s;
          }
          pos >> skipws >> s;              
        }
     Particle p;
   
	bool_point_t periodicityFront = ((Boundary*) m_parent)->periodicityFront();
	MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "periodicityFront = " << periodicityFront);
	bool_point_t periodicityBack = ((Boundary*) m_parent)->periodicityBack();
	MSG_DEBUG("ParticleCreatorWall::adjustBoxSize", "periodicityBack = " 
	<< periodicityBack);
	
        string freeOrFrozen = "free";
        /* peek() looks at what the next character is without extracting is. This loop 
         * continues until the first '!' at the end of the input file is seen. 
         * 
         */
        
        
        while(!pos.eof() && pos.peek() != '!' ){ 
           pos >> skipws >> freeOrFrozen ;  //first species, dismissed.
           readParticle(p,pos,freeOrFrozen); 
           if(!m_particlesInside){        
                if (!(M_BOUNDARY->isInside(p.r))) {
                 
                    for(int i = 0 ; i< SPACE_DIMS; i++){
                        if((p.r[i]-size[i]) <=0){ //smaller than box
                            frameRCfront[i] = frameRCfront[i] || !periodicityFront[i] ;
                        }
                        if((size[i]-p.r[i])<=0){ //bigger than box
                            frameRCend[i] = frameRCend[i] || !periodicityBack[i];
                        }          
                    }
                }
            }   
        
        }
  
        
        	MSG_DEBUG("ParticleCreatorFile::adjustBoxSize","frameRCfront = " << frameRCfront << ", frameRCend = " << frameRCend);

}
