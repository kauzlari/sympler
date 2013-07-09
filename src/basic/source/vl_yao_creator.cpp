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




#include <algorithm>
#include <ctime>

#include "vl_yao_creator.h"
#include "threads.h"
#include "inlet_cell.h"
#include "simulation.h"
#include "vertex_list.h"
#include "manager_cell.h"
#include "particle.h"
#include "particle_cache.h"
#include "particle_creator.h"
#include "colour_pair.h"
#include "phase.h"

#ifdef _OPENMP
  #include "omp.h"
#endif



// necessary for arguments that are pointers to Phase
#include "phase.h"

using namespace std;


#define NUM_NEIGHBORS 26
#define NUM_THIRD_NEIGHBORS 9

#define M_PHASE  ((Phase*) m_parent)
#define M_MANAGER  M_PHASE->manager()
#define M_SIMULATION ((Simulation*) M_PHASE->parent())

const PairCreator_Register<VLYaoCreator> vl_yao_creator("VLYaoCreator");



VLYaoCreator::VLYaoCreator(Phase* p): PairCreator(p)
{
  init();
}


VLYaoCreator::~VLYaoCreator()
{
}

void VLYaoCreator::init()
{
  m_properties.setClassName("VLYaoCreator");

  m_properties.setDescription("Should create fast verlet list pairs but the current implementation is too slow and hence the usage of this module can usually not be recommended. The algorhitm is based on the method, presented by Zhenhua Yao, Jian-Sheng Wang and Min Cheng.");

  DOUBLEPC
    (skinSize,
     m_skin_size,
     0,
     "The size of the additional skin to the cut-off to get the verlet-list cut-off.");

  STRINGPC
    (displacement, m_displacement_name,
     "Name of the displacement.");

  INTPC(slotSize, m_n_slots, 0,
     "The number of slots that are needed for saving the neighbour of a particle.");

  m_displacement_name = "displacement";
  m_create_now = true;
  m_n_slots = 0;
  m_skin_size = 0.;
}


void VLYaoCreator::setup()
{
  PairCreator::setup();

  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

  if (m_n_slots == 0)
    throw gError("VLYaoCreator::setup", "Please define a positive non-zero maximum number of slots for the neighbours of a particle with the attribute 'slotSize'!");

  if (m_skin_size < 0.)
    throw gError("VLYaoCreator::setup", "Please define a positive skin size!");

//   m_cutoff = M_SIMULATION->maxCutoff + m_skin_size;
  FOR_EACH_COLOUR_PAIR
          (M_MANAGER,

	    cp->setCutoff(M_SIMULATION->maxCutoff + m_skin_size);
	  );

  m_capacity_free.resize(M_MANAGER->nColours());
  m_capacity_frozen.resize(M_MANAGER->nColours());
  m_offset_free.resize(M_MANAGER->nColours());
  m_offset_frozen.resize(M_MANAGER->nColours());
  m_offset.resize(M_MANAGER->nColours());
//   vector<size_t> index(M_MANAGER->nColours());

  for (size_t colour = 0; colour < M_MANAGER->nColours(); ++colour) {
    m_capacity_free[colour].resize(M_MANAGER->nColours());
    m_capacity_frozen[colour].resize(M_MANAGER->nColours());
    m_offset_free[colour].resize(M_MANAGER->nColours());
    m_offset_frozen[colour].resize(M_MANAGER->nColours());
  }

  for (size_t colour = 0; colour < M_MANAGER->nColours(); ++colour) {

    if (Particle::s_tag_format[colour].attrExists(m_displacement_name)) {
      DataFormat::attribute_t attrDisp = Particle::s_tag_format[colour].attrByName(m_displacement_name);

      if(attrDisp.datatype != DataFormat::POINT)
        throw gError("VLYaoCreator::setup", "the symbol " + m_displacement_name + " is not registered as a point");

      m_displacement_o = Particle::s_tag_format[colour].attrByName(m_displacement_name).offset;
    }
    else {
      throw gError("VLYaoCreator::setup", "No symbol '" + m_displacement_name + "' found! You need another module introducing it. This might be for example an IntegratorVelocityVerletDisp.");
    }

    string displacementOld_name = string(m_displacement_name + "__Old");

    if (!Particle::s_tag_format[colour].attrExists(displacementOld_name)) {      
      m_displacementOld_o = 
	Particle::s_tag_format[colour].addAttribute
	(displacementOld_name,
	 DataFormat::POINT,
	 true, // must be persistent!
	 displacementOld_name).offset;
    }
    else {
      throw gError("VLYaoCreator::setup", "Symbol '" + displacementOld_name + "' already existing for species '" + M_MANAGER->species(colour) + "' but I have to register mine! You may not define it twice, so please change the other symbol's name!");
    }

    ColourPair *cp = M_MANAGER->cp(colour, colour);
    if (cp->needPairs()) {
      m_capacity_free[colour][colour] = m_n_slots;
      
      if (!ParticleCreator::s_createFrozenParts)
	m_capacity_frozen[colour][colour] = 0;
      else
	m_capacity_frozen[colour][colour] = m_n_slots;
      
      Particle::s_tag_format[colour].addAttribute
        ("size" + ObjToString(colour) + "free",
         DataFormat::INT,
         true,
         "s1");
      
      for (int i = 0; i < m_capacity_free[colour][colour]; i++) {
	Particle::s_tag_format[colour].addAttribute
	  ("slot" + ObjToString(colour) + "free" + ObjToString(i),
	   DataFormat::INT,
	   true,
	   "sl1");
      }
      
      Particle::s_tag_format[colour].addAttribute
        ("size" + ObjToString(colour) + "frozen",
         DataFormat::INT,
         true,
         "s1");
      
      for (int i = 0; i < m_capacity_frozen[colour][colour]; i++) {
	Particle::s_tag_format[colour].addAttribute
	  ("slot" + ObjToString(colour) + "frozen" + ObjToString(i),
	   DataFormat::INT,
	   true,
	   "sl1");
      }
    }
    else {
      m_capacity_free[colour][colour] = 0;
      m_capacity_frozen[colour][colour] = 0;
      Particle::s_tag_format[colour].addAttribute
        ("size" + ObjToString(colour) + "free",
         DataFormat::INT,
         true,
         "s1");
      
      Particle::s_tag_format[colour].addAttribute
        ("size" + ObjToString(colour) + "frozen",
         DataFormat::INT,
         true,
         "s1");
    }
    
    for(size_t c2 = colour+1; c2 < M_MANAGER->nColours(); ++c2) {
      cp = M_MANAGER->cp(colour, c2);
      if (cp->needPairs()) {
	
	m_capacity_free[colour][c2] = m_n_slots;
	m_capacity_frozen[colour][c2] = m_n_slots;
	
	m_capacity_free[c2][colour] = m_n_slots;
	m_capacity_frozen[c2][colour] = m_n_slots;
	
	Particle::s_tag_format[colour].addAttribute
	  ("size" + ObjToString(c2) + "free",
	   DataFormat::INT,
	   true,
	   "s1");
	
	for (int i = 0; i < m_capacity_free[colour][c2]; i++) {
	  Particle::s_tag_format[colour].addAttribute
	    ("slot" + ObjToString(c2) + "free" + ObjToString(i),
	     DataFormat::INT,
	     true,
	     "sl1");
	}
	
	
	Particle::s_tag_format[colour].addAttribute
	  ("size" + ObjToString(c2) + "frozen",
	   DataFormat::INT,
	   true,
	   "s1");
	
	for (int i = 0; i < m_capacity_frozen[colour][c2]; i++) {
	  Particle::s_tag_format[colour].addAttribute
	    ("slot" + ObjToString(c2) + "frozen" + ObjToString(i),
	     DataFormat::INT,
	     true,
	     "sl1");
	}
	
	
	Particle::s_tag_format[c2].addAttribute
	  ("size" + ObjToString(colour) + "free",
	   DataFormat::INT,
	   true,
	   "s1");
	
	for (int i = 0; i < m_capacity_free[c2][colour]; i++) {
	  Particle::s_tag_format[c2].addAttribute
	    ("slot" + ObjToString(colour) + "free" + ObjToString(i),
	     DataFormat::INT,
	     true,
	     "sl1");
	}
	
	Particle::s_tag_format[c2].addAttribute
	  ("size" + ObjToString(colour) + "frozen",
	   DataFormat::INT,
	   true,
	   "s1");
	
	for (int i = 0; i < m_capacity_frozen[c2][colour]; i++) {
	  Particle::s_tag_format[c2].addAttribute
	    ("slot" + ObjToString(colour) + "frozen" + ObjToString(i),
	     DataFormat::INT,
	     true,
	     "sl1");
	}
	
      }
      else {
	m_capacity_free[colour][c2] = 0;
	m_capacity_free[c2][colour] = 0;
	m_capacity_frozen[colour][c2] = 0;
	m_capacity_frozen[c2][colour] = 0;
	
	Particle::s_tag_format[colour].addAttribute
	  ("size" + ObjToString(c2) + "free",
	   DataFormat::INT,
	   true,
	   "s1");
	
	Particle::s_tag_format[colour].addAttribute
	  ("size" + ObjToString(c2) + "frozen",
	   DataFormat::INT,
	   true,
	   "s1");
	
	Particle::s_tag_format[c2].addAttribute
	  ("size" + ObjToString(colour) + "free",
	   DataFormat::INT,
	   true,
	   "s1");
	
	Particle::s_tag_format[c2].addAttribute
	  ("size" + ObjToString(colour) + "frozen",
	   DataFormat::INT,
	   true,
	   "s1");
      }
    }
  }
  
  size_t index = Particle::s_tag_format[0].attrByName("size" + ObjToString(0) + "free").index;

  m_int_offset =
  Particle::s_tag_format[0].attrByIndex(index+1).offset - Particle::s_tag_format[0].attrByIndex(index).offset;

  for (size_t c1 = 0; c1 < M_MANAGER->nColours(); ++c1) {
    m_offset_free[c1][0] = Particle::s_tag_format[c1].attrByName("size" + ObjToString(0) + "free").offset/*m_offset[c1]*/;
    m_offset_frozen[c1][0] = m_offset_free[c1][0] + (m_capacity_free[c1][0] + 1)*m_int_offset;
  }


    for (size_t c1 = 0; c1 < M_MANAGER->nColours(); ++c1) {
      for (size_t c2 = 0; c2 < M_MANAGER->nColours(); ++c2) {
        for (size_t c3 = 0; c3 < c2; ++c3) {
	  m_offset_free[c1][c2] /*+*/= m_offset_free[c1][c3] + ((m_capacity_free[c1][c3] + 1)*m_int_offset + (m_capacity_frozen[c1][c3] + 1)*m_int_offset);
	  m_offset_frozen[c1][c2] /*+*/= m_offset_frozen[c1][c3] + ((m_capacity_free[c1][c3] + 1)*m_int_offset + (m_capacity_frozen[c1][c3] + 1)*m_int_offset); //works this way???
        }
      }
    }

    string m_filename = "output";
    m_s.open(m_filename.c_str());
    m_s.flags(ios::fixed);

    writeHeader();
}

void VLYaoCreator::setupAfterParticleCreation()
{

  // initialise old displacement
  FOR_EACH_PARTICLE
    (M_PHASE,
     for (size_t col = 0; col < M_MANAGER->nColours(); ++col) {
       for(size_t dir = 0; dir < SPACE_DIMS; ++dir)
	 (i->tag.pointByOffset(m_displacementOld_o))[dir] = 0; 
     }
    );

}

void VLYaoCreator::setDtf() {

  //  MSG_DEBUG("VLYaoCreator::setDtf", "start");

  size_t n_colours = M_MANAGER->nColours();

//  The Yao way for setting Cells' dirty flags

  Cell* currentCell;

//   size_t countDirty = 0;

  LL_FOR_EACH__PARALLEL
   (Cell,
    M_MANAGER->firstCell(),
    M_MANAGER->activeCells(),
    NULL,

    if(!i->cellUsed()) {
      double max_disp = 0;
      double max2 = 0;
      
      for (size_t c = 0; c < n_colours; ++c) {
	for (list<Particle*>::iterator j = i->particles(c).begin(); j != i->particles(c).end(); ++j) {
	  point_t dispNow = ((*j)->tag.pointByOffset(this->m_displacement_o))-((*j)->tag.pointByOffset(this->m_displacementOld_o));
	  double tempDisp = dispNow.abs();
	  
	  if (max_disp < tempDisp)
	    max_disp = tempDisp;
	  else if (max2 < tempDisp)
	    max2 = tempDisp;
	}
	
	
	if ((max_disp + max2) > m_skin_size) {
	  m_dtf_cells.push_back(i);
	  i->cellUsed() = true;
//           ++countDirty;
	  // store the current displacements as old displacements for the particles in this Cell
	  for (list<Particle*>::iterator j = i->particles(c).begin(); j != i->particles(c).end(); ++j) {
	    ((*j)->tag.pointByOffset(this->m_displacementOld_o)) = ((*j)->tag.pointByOffset(this->m_displacement_o));
	  }
	  for (int n = 0; n < NUM_NEIGHBORS; ++n) {
	    for (list<CellLink*>::iterator cl = i->neighbors(n).begin(); cl != i->neighbors(n).end(); ++cl) {
	      currentCell = (*cl)->first();
	      if (!currentCell->cellUsed()) {
		m_dtf_cells.push_back(currentCell);
		currentCell->cellUsed() = true;
//                 ++countDirty;
		// store the current displacements as old displacements for the particles in this Cell
		for (list<Particle*>::iterator j = currentCell->particles(c).begin(); j != currentCell->particles(c).end(); ++j) {
		  ((*j)->tag.pointByOffset(this->m_displacementOld_o)) = ((*j)->tag.pointByOffset(this->m_displacement_o));
		}
	      }
	      currentCell = (*cl)->second();
	      if (!currentCell->cellUsed()) {
		m_dtf_cells.push_back(currentCell);
		currentCell->cellUsed() = true;
//                 ++countDirty;
		// store the current displacements as old displacements for the particles in this Cell
		for (list<Particle*>::iterator j = currentCell->particles(c).begin(); j != currentCell->particles(c).end(); ++j) {
		  ((*j)->tag.pointByOffset(this->m_displacementOld_o)) = ((*j)->tag.pointByOffset(this->m_displacement_o));
		}
	      }
	    }    
	  }
	}
      }
    }
    );

//    MSG_DEBUG("VLYaoCreator::setDtf", "# of dirty cells" << countDirty);


// commented out because resetting of displacement was missing 
// and too complicated to implement in this method
#if 0

bool pushCell = false;
FOR_EACH_FREE_PARTICLE
    (M_PHASE,
     point_t disp = i->tag.pointByOffset(this->m_displacement_o);

     if (disp.abs() > 0.5*m_skin_size) {
       pushCell = true;
       break;
     }
    );


#ifdef _OPENMP
  if (pushCell) {
#pragma omp parallel for
    for (int t = 0; t < global::n_threads; ++t) {
      CellLink* first = M_MANAGER->firstLink()[t];	
      for (CellLink* i = first; i != NULL; i = i->next) {   
//	  for (vector<CellLink*>::iterator i = (M_MANAGER->links()[t]).begin(); i != (M_MANAGER->links()[t]).end(); ++i) {

        bool found = false;
        if ((!i->first()->cellUsed()) || (!i->second()->cellUsed())) {
          point_t tmp;
      // If frozen particles exist, we check for single displacements of the other cell
        for (size_t c1 = 0; c1 < n_colours; ++c1) {
          for (size_t c2 = 0; c2 < n_colours; ++c2) {
            if (!(i->first()->frozenParticles(c1).empty())) {
              if (!(i->second()->particles(c2).empty())) {
                for (list<Particle*>::iterator j = i->second()->particles(c2).begin(); j != i->second()->particles(c2).end(); ++j) {
                  if (((*j)->tag.pointByOffset(this->m_displacement_o)).abs() > m_skin_size) {
                    if (!i->first()->cellUsed()) {
                      m_dtf_cells.push_back(i->first());
                      i->first()->cellUsed() = true;
                    }
                    if (!i->second()->cellUsed()) {
                      m_dtf_cells.push_back(i->second());
                      i->second()->cellUsed() = true;
                    }

                    found = true;
		    break;
                  }
                  if (found)
                    break;
              }
            }
          }
            if (!found) {
              if (!(i->second()->frozenParticles(c1).empty())) {
                if (!(i->first()->particles(c2).empty())) {
                  for (list<Particle*>::iterator k = i->first()->particles(c2).begin(); k != i->first()->particles(c2).end(); ++k) {
                    if (((*k)->tag.pointByOffset(this->m_displacement_o)).abs() > m_skin_size) {
                      if (!i->first()->cellUsed()) {
                        m_dtf_cells.push_back(i->first());
                        i->first()->cellUsed() = true;
                      }
                      if (!i->second()->cellUsed()) {
                        m_dtf_cells.push_back(i->second());
                        i->second()->cellUsed() = true;
                      }
                      found = true;
		      break;
                    }
                    if (found)
                      break;
                  }
                }
              }
            }
            if (found)
              break;
          }
          if (found)
            break;
        }

        if (!found) {
          // We check for relative displacements of two cells with free particles
          for (size_t c1 = 0; c1 < n_colours; ++c1) {
            for (size_t c2 = 0; c2 < n_colours; ++c2) {
              for (list<Particle*>::iterator _j = i->first()->particles(c1).begin(); _j != i->first()->particles(c1).end(); ++_j) {
                point_t dispOne = ((*_j)->tag.pointByOffset(this->m_displacement_o));

                for (list<Particle*>::iterator _k = i->second()->particles(c2).begin(); _k != i->second()->particles(c2).end(); ++_k) {
                  point_t dispTwo = ((*_k)->tag.pointByOffset(this->m_displacement_o));

                  for (int s = 0; s < SPACE_DIMS; s++)
                    tmp[s] = (dispOne[s] - dispTwo[s]);

                  if (tmp.abs() > m_skin_size) {
                    if (!i->first()->cellUsed()) {
                      m_dtf_cells.push_back(i->first());
                      i->first()->cellUsed() = true;
                    }
                    if (!i->second()->cellUsed()) {
                      m_dtf_cells.push_back(i->second());
                      i->second()->cellUsed() = true;
                    }

                    found = true;
                    break;
                  }
                }
                if (found)
                break;
              }
              if (found)
                break;
            } // Loop over c2
            if (found)
              break;
          } // Loop over c1
        }

      }
    }
  }
}

#else
  if (pushCell) {
    LL_FOR_EACH__PARALLEL
    (CellLink,
      M_MANAGER->firstLink(),
      M_MANAGER->activeLinks(),
      NULL,

    bool found = false;
    if ((!i->first()->cellUsed()) || (!i->second()->cellUsed())) {
      point_t tmp;
      // If frozen particles exist, we check for single displacements of the other cell
      for (size_t c1 = 0; c1 < n_colours; ++c1) {
        for (size_t c2 = 0; c2 < n_colours; ++c2) {
          if (!(i->first()->frozenParticles(c1).empty())) {
            if (!(i->second()->particles(c2).empty())) {
              for (list<Particle*>::iterator j = i->second()->particles(c2).begin(); j != i->second()->particles(c2).end(); ++j) {
                if (((*j)->tag.pointByOffset(this->m_displacement_o)).abs() > m_skin_size) {
                  if (!i->first()->cellUsed()) {
                    m_dtf_cells.push_back(i->first());
                    i->first()->cellUsed() = true;
                  }
                  if (!i->second()->cellUsed()) {
                    m_dtf_cells.push_back(i->second());
                    i->second()->cellUsed() = true;
                  }

                  found = true;
		  break;
                }
                if (found)
                  break;
            }
          }
        }
          if (!found) {
            if (!(i->second()->frozenParticles(c1).empty())) {
              if (!(i->first()->particles(c2).empty())) {
                for (list<Particle*>::iterator k = i->first()->particles(c2).begin(); k != i->first()->particles(c2).end(); ++k) {
                  if (((*k)->tag.pointByOffset(this->m_displacement_o)).abs() > m_skin_size) {
                    if (!i->first()->cellUsed()) {
                      m_dtf_cells.push_back(i->first());
                      i->first()->cellUsed() = true;
                    }
                    if (!i->second()->cellUsed()) {
                      m_dtf_cells.push_back(i->second());
                      i->second()->cellUsed() = true;
                    }
                    found = true;
		    break;
                  }
                  if (found)
                    break;
                }
              }
            }
          }
          if (found)
            break;
        }
        if (found)
        break;
      }

      if (!found) {
        // We check for relative displacements of two cells with free particles
        for (size_t c1 = 0; c1 < n_colours; ++c1) {
          for (size_t c2 = 0; c2 < n_colours; ++c2) {
            for (list<Particle*>::iterator _j = i->first()->particles(c1).begin(); _j != i->first()->particles(c1).end(); ++_j) {
              point_t dispOne = ((*_j)->tag.pointByOffset(this->m_displacement_o));

              for (list<Particle*>::iterator _k = i->second()->particles(c2).begin(); _k != i->second()->particles(c2).end(); ++_k) {
                point_t dispTwo = ((*_k)->tag.pointByOffset(this->m_displacement_o));

                for (int s = 0; s < SPACE_DIMS; s++)
                  tmp[s] = (dispOne[s] - dispTwo[s]);

                if (tmp.abs() > m_skin_size) {
                  if (!i->first()->cellUsed()) {
                    m_dtf_cells.push_back(i->first());
                    i->first()->cellUsed() = true;
                  }
                  if (!i->second()->cellUsed()) {
                    m_dtf_cells.push_back(i->second());
                    i->second()->cellUsed() = true;
                  }

                  found = true;
                  break;
                }
              }
              if (found)
                break;
            }
            if (found)
              break;
          } // Loop over c2
          if (found)
            break;
        } // Loop over c1
      }

    }
    );
  }
#endif
  
#endif // end of commented out method

//   MSG_DEBUG("VLYaoCreator::setDtf", "done");

}


// Here not all directions have to be checked because in the same cell either everything is deleted or nothing.
void VLYaoCreator::pairAddedFree(Pairdist* pd)
{
  if (pd->firstPart()->r[0] < pd->secondPart()->r[0]) {
    int* size = &(pd->firstPart()->tag.intByOffset(m_offset_free[pd->firstPart()->c][pd->secondPart()->c]));
    (*size) +=1;
//     pd->firstPart()->tag.intByOffset(m_offset_free[pd->firstPart()->c][pd->secondPart()->c]) += 1;

    if ((*size) < m_capacity_free[pd->firstPart()->c][pd->secondPart()->c]) {
      pd->firstPart()->tag.intByOffset(m_offset_free[pd->firstPart()->c][pd->secondPart()->c] + (*size)*m_int_offset) = pd->mySlot;
    }
    else {
      throw gError("VLYaoCreator::pairAddedFree", "There is not enough space to save all pairs of the current type. Please define a bigger slotSize! Reached number of pairs: " + ObjToString(*size) + ", capacity = " + ObjToString(m_capacity_free[pd->firstPart()->c][pd->secondPart()->c]));
    }
  }

  //Add the pair to the tag of the second Particle else
  else {
    int* size = &(pd->secondPart()->tag.intByOffset(m_offset_free[pd->secondPart()->c][pd->firstPart()->c]));
    (*size) +=1;
//     pd->secondPart()->tag.intByOffset(m_offset_free[pd->secondPart()->c][pd->firstPart()->c]) += 1;

    if ((*size) < m_capacity_free[pd->secondPart()->c][pd->firstPart()->c]) {
      pd->secondPart()->tag.intByOffset(m_offset_free[pd->secondPart()->c][pd->firstPart()->c] + (*size)*m_int_offset) = pd->mySlot;
    }
    else {
      throw gError("VLYaoCreator::pairAddedFree", "There is not enough space to save all pairs of the current type. Please define a bigger slotSize! Reached number of pairs: " + ObjToString(*size) + ", capacity = " + ObjToString(m_capacity_free[pd->secondPart()->c][pd->firstPart()->c]));
    }

  }
}


void VLYaoCreator::pairAddedFrozen(Pairdist* pd)
{
  //Add the pair to the tag of the first Particle if the x-position of the first is smaller than the x-position of the second
  if (pd->firstPart()->r[0] < pd->secondPart()->r[0]) {
    pd->firstPart()->tag.intByOffset(m_offset_frozen[pd->firstPart()->c][pd->secondPart()->c]) += 1;
    int size = pd->firstPart()->tag.intByOffset(m_offset_frozen[pd->firstPart()->c][pd->secondPart()->c]);

      if (size < m_capacity_frozen[pd->firstPart()->c][pd->secondPart()->c]) {
	pd->firstPart()->tag.intByOffset(m_offset_frozen[pd->firstPart()->c][pd->secondPart()->c] + size*m_int_offset) = pd->mySlot;
      }
      else {
	throw gError("VLYaoCreator::pairAddedFrozen", "There is not enough space to save all pairs of the current type. Please define a bigger slotSize!");
      }
  }

  //Add the pair to the tag of the second Particle else
  else {
    pd->secondPart()->tag.intByOffset(m_offset_frozen[pd->secondPart()->c][pd->firstPart()->c]) += 1;
    int size = pd->secondPart()->tag.intByOffset(m_offset_frozen[pd->secondPart()->c][pd->firstPart()->c]);

      if (size < m_capacity_frozen[pd->secondPart()->c][pd->firstPart()->c]) {
        pd->secondPart()->tag.intByOffset(m_offset_frozen[pd->secondPart()->c][pd->firstPart()->c] + size*m_int_offset) = pd->mySlot;
      }
      else {
	throw gError("VLYaoCreator::pairAddedFrozen", "There is not enough space to save all pairs of the current type. Please define a bigger slotSize!");
      }
  }
}


void VLYaoCreator::invalidatePositions()
{
  // avoid pair-computation if there are no non-bonded pairs needed
  m_valid_dist = true;

  FOR_EACH_COLOUR_PAIR
  (M_MANAGER,
   m_valid_dist = !(cp->needPairs()) && m_valid_dist;
  );

}


void VLYaoCreator::createDistances1stDirty
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist) {
      list<Particle*>::iterator p1_end = first_p.end();
      list<Particle*>::iterator p2_end = second_p.end();
      for (list<Particle*>::iterator i = first_p.begin(); i != p1_end; ++i) {
        for (list<Particle*>::iterator j = second_p.begin(); j != p2_end; ++j) {

          dist_t d;

          for (int _i = 0; _i < SPACE_DIMS; _i++) {
            d.cartesian[_i] = -dir * cell_dist[_i]
              + (*i) -> r[_i] - first_c -> corner1[_i]
              - (*j) -> r[_i] + second_c -> corner1[_i];
          }




          bool create = false;
//           if ((*i)->r.x < (*j)->r.x) {
          if (d.cartesian.x < 0) {
            create = true;

          }
//           else if ((*i)->r.x == (*j)->r.x) {
          else if (d.cartesian.x == 0) {
//             if ((*i)->r.y < (*j)->r.y) {
            if (d.cartesian.y < 0) {
              create = true;

            }
//             else if ((*i)->r.y == (*j)->r.y) {
            else if (d.cartesian.y == 0) {
//               if ((*i)->r.z < (*j)->r.z) {
              if (d.cartesian.z < 0) {
                create = true;

              }
//               else if ((*i)->r.z == (*j)->r.z)
              else if (d.cartesian.z == 0)
                throw gError("VLYaoCreator::createDistances1stDirty", " Positions of particle " + ObjToString((*i)->mySlot) + " and " + ObjToString((*j)->mySlot) + " are identical. ");
            }
          }
          if (create) {
            d.abs_square = 0;

            for (int _i = 0; _i < SPACE_DIMS; _i++) {
              d.abs_square += d.cartesian[_i]*d.cartesian[_i];
            }
            /* Take care: The order of *j, *i defines the direction d \
              is pointing to. */
            if (d.abs_square < cutoff_sq) {
              d.abs = sqrt(d.abs_square);
              Pairdist* temp = &distances[PairCreator::counterTN].newPair();
              temp->set(d, *i, *j, ao_f, ao_s);

              (*i)->tag.intByOffset(m_offset_free[(*i)->c][(*j)->c]) += 1;
              int _size = (*i)->tag.intByOffset(m_offset_free[(*i)->c][(*j)->c]);
              if (_size < m_capacity_free[(*i)->c][(*j)->c]) {
                (*i)->tag.intByOffset(m_offset_free[(*i)->c][(*j)->c] + _size*m_int_offset) = temp->mySlot;
              }
              else {
                throw gError("VLYaoCreator::createDistances1stDirty", "There is not enough space to save all pairs of colour " + ObjToString((*j)->c) + " for particle " + ObjToString((*i)->mySlot) + " of colour " + ObjToString((*i)->c) + ". Please define a bigger slotSize!");
              }

              ++PairCreator::counterTN;
              if (PairCreator::counterTN == global::n_threads)
                PairCreator::counterTN = 0;
            }
          }
        }
      }
}


void VLYaoCreator::createDistances1stDirtyFrozen
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist) {
      list<Particle*>::iterator p1_end = first_p.end();
      list<Particle*>::iterator p2_end = second_p.end();
      for (list<Particle*>::iterator i = first_p.begin(); i != p1_end; ++i) {
        for (list<Particle*>::iterator j = second_p.begin(); j != p2_end; ++j) {

          dist_t d;

          for (int _i = 0; _i < SPACE_DIMS; _i++) {
            d.cartesian[_i] = -dir * cell_dist[_i]
              + (*i) -> r[_i] - first_c -> corner1[_i]
              - (*j) -> r[_i] + second_c -> corner1[_i];
          }




          bool create = false;
//           if ((*i)->r.x < (*j)->r.x) {
          if (d.cartesian.x < 0) {
            create = true;

          }
//           else if ((*i)->r.x == (*j)->r.x) {
          else if (d.cartesian.x == 0) {
//             if ((*i)->r.y < (*j)->r.y) {
            if (d.cartesian.y < 0) {
              create = true;

            }
//             else if ((*i)->r.y == (*j)->r.y) {
            else if (d.cartesian.y == 0) {
//               if ((*i)->r.z < (*j)->r.z) {
              if (d.cartesian.z < 0) {
                create = true;

              }
//               else if ((*i)->r.z == (*j)->r.z)
              else if (d.cartesian.z == 0)
                throw gError("VLYaoCreator::createDistances1stDirty", " Positions of particle " + ObjToString((*i)->mySlot) + " and " + ObjToString((*j)->mySlot) + " are identical. ");
            }
          }
          if (create) {
            d.abs_square = 0;

            for (int _i = 0; _i < SPACE_DIMS; _i++) {
              d.abs_square += d.cartesian[_i]*d.cartesian[_i];
            }
            /* Take care: The order of *j, *i defines the direction d \
              is pointing to. */
            if (d.abs_square < cutoff_sq) {
              d.abs = sqrt(d.abs_square);
              Pairdist* temp = &distances[PairCreator::counterTN].newPair();
              temp->set(d, *i, *j, ao_f, ao_s);

              (*i)->tag.intByOffset(m_offset_frozen[(*i)->c][(*j)->c]) += 1;
              int _size = (*i)->tag.intByOffset(m_offset_frozen[(*i)->c][(*j)->c]);
              if (_size < m_capacity_frozen[(*i)->c][(*j)->c]) {
                (*i)->tag.intByOffset(m_offset_frozen[(*i)->c][(*j)->c] + _size*m_int_offset) = temp->mySlot;
              }
              else {
                throw gError("VLYaoCreator::createDistances1stDirtyFrozen", "There is not enough space to save all pairs of colour " + ObjToString((*j)->c) + " for particle " + ObjToString((*i)->mySlot) + " of colour " + ObjToString((*i)->c) + ". Please define a bigger slotSize!");
              }

              ++PairCreator::counterTN;
              if (PairCreator::counterTN == global::n_threads)
                PairCreator::counterTN = 0;
            }
          }
        }
      }
}


void VLYaoCreator::createDistances2ndDirty
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist) {
      list<Particle*>::iterator p2_end = second_p.end();
      list<Particle*>::iterator p1_end = first_p.end();
      for (list<Particle*>::iterator i = second_p.begin(); i != p2_end; ++i) {
        for (list<Particle*>::iterator j = first_p.begin(); j != p1_end; ++j) {

          dist_t d;

          for (int _i = 0; _i < SPACE_DIMS; _i++) {
            d.cartesian[_i] = -dir * cell_dist[_i]
              + (*j) -> r[_i] - first_c -> corner1[_i]
              - (*i) -> r[_i] + second_c -> corner1[_i];
          }
// if (distance == frozenPairs(0))
// MSG_DEBUG("VLYaoCreator::createDistances2ndDirty", " slot of j = " << (*j)->mySlot);
          bool create = false;
//           if ((*i)->r.x < (*j)->r.x) {
          if (d.cartesian.x > 0) {
            create = true;

          }
//           else if ((*i)->r.x == (*j)->r.x) {
          else if (d.cartesian.x == 0) {
//             if ((*i)->r.y < (*j)->r.y) {
            if (d.cartesian.y > 0) {
              create = true;

            }
//             else if ((*i)->r.y == (*j)->r.y) {
            else if (d.cartesian.y == 0) {
//               if ((*i)->r.z < (*j)->r.z) {
              if (d.cartesian.z > 0) {
                create = true;

              }
//               else if ((*i)->r.z == (*j)->r.z)
              else if (d.cartesian.z == 0)
                throw gError("VLYaoCreator::createDistances2ndDirty", " Positions of particle " + ObjToString((*i)->mySlot) + " and " + ObjToString((*j)->mySlot) + " are identical. ");
            }
          }
          if (create) {
            d.abs_square = 0;

            for (int _i = 0; _i < SPACE_DIMS; _i++) {
              d.abs_square += d.cartesian[_i]*d.cartesian[_i];
            }

            /* Take care: The order of *j, *i defines the direction d \
              is pointing to. */
            if (d.abs_square < cutoff_sq) {
              d.abs = sqrt(d.abs_square);
              Pairdist* temp = &distances[PairCreator::counterTN].newPair();
              temp->set(d, *j, *i, ao_f, ao_s);

              (*i)->tag.intByOffset(m_offset_free[(*i)->c][(*j)->c]) += 1;
              int _size = (*i)->tag.intByOffset(m_offset_free[(*i)->c][(*j)->c]);
              if (_size < m_capacity_free[(*i)->c][(*j)->c]) {
                (*i)->tag.intByOffset(m_offset_free[(*i)->c][(*j)->c] + _size*m_int_offset) = temp->mySlot;
              }
              else {
                throw gError("VLYaoCreator::createDistances1stDirty", "There is not enough space to save all pairs of colour " + ObjToString((*j)->c) + " for particle " + ObjToString((*i)->mySlot) + " of colour " + ObjToString((*i)->c) + ". Please define a bigger slotSize!");
              }

              ++PairCreator::counterTN;
              if (PairCreator::counterTN == global::n_threads)
                PairCreator::counterTN = 0;
            }
          }
        }
      }
}


void VLYaoCreator::createDistances2ndDirtyFrozen
  (vector<PairList> &distances,
   double cutoff_sq,
   int dir,
   Cell *first_c,
   Cell *second_c,
   list<Particle*> &first_p,
   list<Particle*> &second_p,
   bool ao_f,
   bool ao_s,
   point_t &cell_dist) {
      list<Particle*>::iterator p2_end = second_p.end();
      list<Particle*>::iterator p1_end = first_p.end();
      for (list<Particle*>::iterator i = second_p.begin(); i != p2_end; ++i) {
        for (list<Particle*>::iterator j = first_p.begin(); j != p1_end; ++j) {

// if(/*(*j)->mySlot == 62 && */(*i)->mySlot == 46 || (*j)->mySlot == 46)
// MSG_DEBUG("", " pos of part i = " << (*i)->r << " pos of j = " << (*j)->r << " slots for i = " << (*i)->mySlot << " slot of j = " << (*j)->mySlot);
          dist_t d;
          for (int _i = 0; _i < SPACE_DIMS; _i++) {
            d.cartesian[_i] = -dir * cell_dist[_i]
              + (*j) -> r[_i] - first_c -> corner1[_i]
              - (*i) -> r[_i] + second_c -> corner1[_i];
          }
// if (distance == frozenPairs(0))
// MSG_DEBUG("VLYaoCreator::createDistances2ndDirty", " slot of j = " << (*j)->mySlot);
          bool create = false;
//           if ((*i)->r.x < (*j)->r.x) {
          if (d.cartesian.x > 0) {
            create = true;
          }
//           else if ((*i)->r.x == (*j)->r.x) {
          else if (d.cartesian.x == 0) {
//             if ((*i)->r.y < (*j)->r.y) {
            if (d.cartesian.y > 0) {
              create = true;

            }
//             else if ((*i)->r.y == (*j)->r.y) {
            else if (d.cartesian.y == 0) {
//               if ((*i)->r.z < (*j)->r.z) {
              if (d.cartesian.z > 0) {
                create = true;

              }
//               else if ((*i)->r.z == (*j)->r.z)
              else if (d.cartesian.z == 0)
                throw gError("VLYaoCreator::createDistances2ndDirty", " Positions of particle " + ObjToString((*i)->mySlot) + " and " + ObjToString((*j)->mySlot) + " are identical. ");
            }
          }
          if (create) {
            d.abs_square = 0;

            for (int _i = 0; _i < SPACE_DIMS; _i++) {
              d.abs_square += d.cartesian[_i]*d.cartesian[_i];
            }

            /* Take care: The order of *j, *i defines the direction d \
              is pointing to. */
            if (d.abs_square < cutoff_sq) {
              d.abs = sqrt(d.abs_square);
              Pairdist* temp = &distances[PairCreator::counterTN].newPair();
              temp->set(d, *j, *i, ao_f/*s*/, ao_s/*f*/);

              (*i)->tag.intByOffset(m_offset_frozen[(*i)->c][(*j)->c]) += 1;
              int _size = (*i)->tag.intByOffset(m_offset_frozen[(*i)->c][(*j)->c]);
              if (_size < m_capacity_frozen[(*i)->c][(*j)->c]) {
                (*i)->tag.intByOffset(m_offset_frozen[(*i)->c][(*j)->c] + _size*m_int_offset) = temp->mySlot;
              }
              else {
                throw gError("VLYaoCreator::createDistances1stDirtyFrozen", "There is not enough space to save all pairs of colour " + ObjToString((*j)->c) + " for particle " + ObjToString((*i)->mySlot) + " of colour " + ObjToString((*i)->c) + ". Please define a bigger slotSize!");
              }

              ++PairCreator::counterTN;
              if (PairCreator::counterTN == global::n_threads)
                PairCreator::counterTN = 0;
            }
          }
        }
      }
}


void VLYaoCreator::createDistances()
{
  int thread_no = 0;
  size_t n_colours = M_MANAGER->nColours();
  
  time_t t0, t1, t2, t3, t4, t6, t7, t8, t9;
  time(&t0);
  time(&t8);
  
  if (!m_valid_dist) {
    
    if (M_PHASE->particlesAssigned() == true) {
      LL_FOR_EACH__PARALLEL
	(Cell,
	 M_MANAGER->firstCell(),
	 M_MANAGER->activeCells(),
	 NULL,
	 
	 if (m_create_now) {
	   m_dtf_cells.push_back(i);
	 }
	 else {
	   i->cellUsed() = false;
	 }
	 );
      
      if (m_create_now) {
	
	FOR_EACH_PARTICLE
	  (M_PHASE,
	   for (size_t col = 0; col < M_MANAGER->nColours(); ++col) {
	     // MSG_DEBUG("VLYaoCreator::createDistances", " i->c = " << i->c << " col = " << col);
	     i->tag.intByOffset(m_offset_free[i->c][col]) = 0;
	     i->tag.intByOffset(m_offset_free[col][i->c]) = 0;
	     i->tag.intByOffset(m_offset_frozen[i->c][col]) = 0;
	     i->tag.intByOffset(m_offset_frozen[col][i->c]) = 0;
	   }
	   );
	
      }
      
      if (!m_create_now) {
	m_dtf_cells.clear();
	setDtf();
	
	// time(&t9);
	// MSG_INFO("VLYaoCreator", "SetDTF took  " << difftime(t8, t9) << " seconds");
	// time(&t2);
	
	// Switch free and frozen pairs saved in free particles
	FOR_EACH_FREE_PARTICLE
	  (M_PHASE,
	   for (size_t c1 = 0; c1 < n_colours; ++c1) {
	     ColourPair* cp = M_MANAGER->cp(i->c, c1);
	     
	     for (int thread_no = 0; thread_no < (int)global::n_threads; ++thread_no) {
	       updateAndSwitch(i, i->tag.intByOffset(m_offset_free[i->c][c1])/*size*//*, p->tag.intByOffset(m_offset_free[c1][i->c])*//*sizeP*/, m_offset_free[i->c][c1] /*size_loc*/, m_offset_free[c1][i->c] /*size_locP*/, m_capacity_free[c1][i->c], &cp->freePairs()[thread_no]);
	     }
	   }
	   );
	
	for (size_t c1 = 0; c1 < n_colours; ++c1) {
	  for (size_t c2 = 0; c2 < n_colours; ++c2) {
	    ColourPair* cp = M_MANAGER->cp(c1, c2);
	    
	    for (int thread_no = 0; thread_no < (int)global::n_threads; ++thread_no) {
	      if (!cp->frozenPairs()[thread_no].size() == 0) {
		FOR_EACH_FREE_PARTICLE_C
		  (M_PHASE,
		   c1,
		   
		   updateAndSwitch(i, i->tag.intByOffset(m_offset_frozen[c1][c2])/*size*//*, p->tag.intByOffset(m_offset_frozen[c1][i->c])*//*sizeP*/, m_offset_frozen[c1][c2] /*size_loc*/, m_offset_frozen[c2][c1] /*size_locP*/, m_capacity_frozen[c2][c1], &cp->frozenPairs()[thread_no]);
		   );
		
		FOR_EACH_FROZEN_PARTICLE
		  (M_PHASE,
		   c1,
		   
		   updateAndSwitch(i, i->tag.intByOffset(m_offset_frozen[c1][c2])/*size*//*, p->tag.intByOffset(m_offset_frozen[c1][i->c])*//*sizeP*/, m_offset_frozen[c1][c2] /*size_loc*/, m_offset_frozen[c2][c1] /*size_locP*/, m_capacity_frozen[c2][c1], &cp->frozenPairs()[thread_no]);
		   );
	      }
	    }
	  }
	}
	
	// time(&t3);
	// MSG_INFO("VLYaoCreator", "Switching and updating took " << difftime(t2, t3) << " seconds");
	
	
	// Deleting pairs
	for (list<Cell*>::iterator i = m_dtf_cells.begin(); i != m_dtf_cells.end(); ++i) {
	  size_t n_col = M_MANAGER->nColours();
	  for (size_t _c1 = 0; _c1 < n_col; ++_c1) {
	    for (list<Particle*>::iterator p = (*i)->particles(_c1).begin(); p != (*i)->particles(_c1).end(); ++p) {
	      for (size_t _c2 = 0; _c2 < n_col; ++_c2) {
		ColourPair* cp = M_MANAGER->cp(_c1, _c2);
		
		// Free pairs
		if ((*p)->tag.intByOffset(m_offset_free[_c1][_c2]) > 0) {
		  int size = (*p)->tag.intByOffset(m_offset_free[_c1][_c2]);
		  for (int i1 = 1; i1 <= size; ++i1) {
		    int pair2del = (*p)->tag.intByOffset(m_offset_free[_c1][_c2] + i1*m_int_offset);
		    cp->freePairs()[thread_no].deleteEntry(pair2del);
		  }
		  (*p)->tag.intByOffset(m_offset_free[_c1][_c2]) = 0;
		}
		
		// Frozen pairs
		if ((*p)->tag.intByOffset(m_offset_frozen[_c1][_c2]) > 0) {
		  int size = (*p)->tag.intByOffset(m_offset_frozen[_c1][_c2]);
		  for (int _i1 = 1; _i1 < size+1; ++_i1) {
		    int pair2del = (*p)->tag.intByOffset(m_offset_frozen[_c1][_c2] + _i1*m_int_offset);
		    cp->frozenPairs()[thread_no].deleteEntry(pair2del);
		  }
		  (*p)->tag.intByOffset(m_offset_frozen[_c1][_c2]) = 0;
		}
	      }
	    }
	    
	    for (list<Particle*>::iterator _p = (*i)->frozenParticles(_c1).begin(); _p != (*i)->frozenParticles(_c1).end(); ++_p) {
	      for (size_t __c2 = 0; __c2 < n_col; ++__c2) {
		ColourPair* _cp = M_MANAGER->cp(_c1, __c2);
		
		// Free pairs
		if ((*_p)->tag.intByOffset(m_offset_free[_c1][__c2]) > 0) {
		  int size = (*_p)->tag.intByOffset(m_offset_free[_c1][__c2]);
		  for (int i1 = 1; i1 <= size; ++i1) {
		    int pair2del = (*_p)->tag.intByOffset(m_offset_free[_c1][__c2] + i1*m_int_offset);
		    _cp->freePairs()[thread_no].deleteEntry(pair2del);
		  }
		  (*_p)->tag.intByOffset(m_offset_free[_c1][__c2]) = 0;
		}
		
		// Frozen pairs
		if ((*_p)->tag.intByOffset(m_offset_frozen[_c1][__c2]) > 0) {
		  int size = (*_p)->tag.intByOffset(m_offset_frozen[_c1][__c2]);
		  for (int _i1 = 1; _i1 < size+1; ++_i1) {
		    int pair2del = (*_p)->tag.intByOffset(m_offset_frozen[_c1][__c2] + _i1*m_int_offset);
		    _cp->frozenPairs()[thread_no].deleteEntry(pair2del);
		  }
		  (*_p)->tag.intByOffset(m_offset_frozen[_c1][__c2]) = 0;
		}
	      }
	    }
	  }
	}
      }
      
      time(&t6);
      
      //Creating new pairs
      for (list<Cell*>::iterator dtCell = m_dtf_cells.begin(); dtCell != m_dtf_cells.end(); ++dtCell) {
	/* Loop over the first colour */
        for (size_t c1 = 0; c1 < n_colours; ++c1) {
          /* --- Pairs of the same colour -------------------------------------------------------------- */
          ColourPair *cp = M_MANAGER->cp(c1, c1);
	  
	  // MSG_DEBUG("VLYaoCreator::create pairs", " c = " << c1 << " pointer to cell = " << *dtCell);
	  // MSG_DEBUG("VLYaoCreator::create pairs", " size of parts in DT cell = " << (*dtCell)->particles(c1).size());
          list<Particle*> &dirtyP = (*dtCell)->particles(c1);
          list<Particle*> &dirtyPFrozen = (*dtCell)->frozenParticles(c1);
	  
          if (cp->needPairs()) {
            double cutoff_sq = cp->cutoff() * cp->cutoff();
	    
	    //           list<Particle*>::iterator p_end = dirtyP.end();
	    
            // createDistances for free - free
            for (list<Particle*>::iterator i = dirtyP.begin(); i != dirtyP.end(); ++i) {
              list<Particle*>::iterator j = i;
	      
              for (++j; j != dirtyP.end(); ++j) {
                dist_t d;
                d.abs_square = 0;
                for (int _i = 0; _i < SPACE_DIMS; _i++) {
                  d.cartesian[_i] = (*i)->r[_i] - (*j)->r[_i];
                  d.abs_square += d.cartesian[_i]*d.cartesian[_i];
                }
                /* Take care: The order of *j, *i defines the direction d \
                  is pointing to. */
                if (d.abs_square < cutoff_sq) {
                  d.abs = sqrt(d.abs_square);
                  Pairdist* temp = &cp->freePairs()[PairCreator::counterTN].newPair();
                  temp->set(d, *i, *j, true, true);
                  pairAddedFree(temp);

                  ++PairCreator::counterTN;
                  if (PairCreator::counterTN == global::n_threads)
                    PairCreator::counterTN = 0;
                }
              }
            } // END: createDistances for free - free

            // createDistances for free - frozen
  //           list<Particle*>::iterator p1_end = dirtyP.end();
  //           list<Particle*>::iterator p2_end = dirtyPFrozen.end();

            for (list<Particle*>::iterator i = dirtyP.begin(); i != dirtyP.end(); ++i) {
              for (list<Particle*>::iterator j = dirtyPFrozen.begin(); j != dirtyPFrozen.end(); ++j) {
                dist_t d;
                d.abs_square = 0;
                for (int _i = 0; _i < SPACE_DIMS; _i++) {
                  d.cartesian[_i] = (*i)->r[_i] - (*j)->r[_i];
                  d.abs_square += d.cartesian[_i]*d.cartesian[_i];
                  }

                /* Take care: The order of *j, *i defines the direction d \
                  is pointing to. */
                if (d.abs_square < cutoff_sq) {
                  d.abs = sqrt(d.abs_square);
                  Pairdist* temp = &cp->frozenPairs()[PairCreator::counterTN].newPair();
                  temp->set(d, *i, *j, true, false);

                  pairAddedFrozen(temp);

                  ++PairCreator::counterTN;
                  if (PairCreator::counterTN == global::n_threads)
                    PairCreator::counterTN = 0;
                }
              }
            } // END: createDistances for free - frozen
          }

        for (size_t c2 = c1+1; c2 < n_colours; ++c2) {
          /* --- Pairs of different colour -------------------------------------------------------------- */
          cp = M_MANAGER->cp(c1, c2);
  //         list<Particle*> &dirtyP1 = (*dtCell)->particles(c1);
          list<Particle*> &dirtyP2 = (*dtCell)->particles(c2);
          list<Particle*> &dirtyPFrozen1 = (*dtCell)->frozenParticles(c2);
  //         list<Particle*> &dirtyPFrozen2 = (*dtCell)->frozenParticles(c1);
            // createDistances free(c1) - free(c2)
            if (cp->needPairs()) {
              double cutoff_sq = cp->cutoff() * cp->cutoff();
  //             list<Particle*>::iterator p1_end = dirtyP/*1*/.end();
  //             list<Particle*>::iterator p2_end = dirtyP2.end();
  //             list<Particle*>::iterator p3_end = dirtyPFrozen1.end();
  //             list<Particle*>::iterator p4_end = dirtyPFrozen.end();

              for (list<Particle*>::iterator i = dirtyP/*1*/.begin(); i != dirtyP.end(); ++i) {
                for (list<Particle*>::iterator j = dirtyP2.begin(); j != dirtyP2.end(); ++j) {
                  dist_t d;
                  d.abs_square = 0;
                  for (int _i = 0; _i < SPACE_DIMS; _i++) {
                    d.cartesian[_i] = (*i)->r[_i] - (*j)->r[_i];
                    d.abs_square += d.cartesian[_i]*d.cartesian[_i];
                    }

                  /* Take care: The order of *j, *i defines the direction d \
                    is pointing to. */
                  if (d.abs_square < cutoff_sq) {
                    d.abs = sqrt(d.abs_square);
                    Pairdist* temp = &cp->freePairs()[PairCreator::counterTN].newPair();
                    temp->set(d, *i, *j, true, true);

                    pairAddedFree(temp);

                    ++PairCreator::counterTN;
                    if (PairCreator::counterTN == global::n_threads)
                      PairCreator::counterTN = 0;
                  }
                }
              } // END: createDistances free(c1) - free(c2)

            // createDistances free(c1) - frozen(c2)
            for (list<Particle*>::iterator i = dirtyP/*1*/.begin(); i != dirtyP.end(); ++i) {
              for (list<Particle*>::iterator j = dirtyPFrozen1.begin(); j != dirtyPFrozen1.end(); ++j) {
                dist_t d;
                d.abs_square = 0;
                for (int _i = 0; _i < SPACE_DIMS; _i++) {
                  d.cartesian[_i] = (*i)->r[_i] - (*j)->r[_i];
                  d.abs_square += d.cartesian[_i]*d.cartesian[_i];
                  }
                /* Take care: The order of *j, *i defines the direction d \
                  is pointing to. */
                if (d.abs_square < cutoff_sq) {
                  d.abs = sqrt(d.abs_square);
                  Pairdist* temp = &cp->frozenPairs()[PairCreator::counterTN].newPair();
                  temp->set(d, *i, *j, true, false);

                  pairAddedFrozen(temp);

                  ++PairCreator::counterTN;
                  if (PairCreator::counterTN == global::n_threads)
                    PairCreator::counterTN = 0;
                }
              }
            } // END: createDistances free(c1) - frozen(c2)

            // createDistances frozen(c1) - free(c2)
            for (list<Particle*>::iterator i = dirtyPFrozen/*2*/.begin(); i != dirtyPFrozen.end(); ++i) {
              for (list<Particle*>::iterator j = dirtyP2.begin(); j != dirtyP2.end(); ++j) {
                dist_t d;
                d.abs_square = 0;
                for (int _i = 0; _i < SPACE_DIMS; _i++) {
                  d.cartesian[_i] = (*i)->r[_i] - (*j)->r[_i];
                  d.abs_square += d.cartesian[_i]*d.cartesian[_i];
                  }

                /* Take care: The order of *j, *i defines the direction d \
                  is pointing to. */
                if (d.abs_square < cutoff_sq) {
                  d.abs = sqrt(d.abs_square);
                  Pairdist* temp = &cp->frozenPairs()[PairCreator::counterTN].newPair();
                  temp->set(d, *i, *j, false, true);

                  pairAddedFrozen(temp);

                  ++PairCreator::counterTN;
                  if (PairCreator::counterTN == global::n_threads)
                    PairCreator::counterTN = 0;
                }
              }
            } // END: createDistances frozen(c1) - free(c2)
          }
        } /* Loop over c2 */
      } /* Loop over c1 */

    for (int n = NUM_THIRD_NEIGHBORS; n < NUM_NEIGHBORS; ++n) {
      for (list<CellLink*>::iterator cl = (*dtCell)->neighbors(n).begin(); cl != (*dtCell)->neighbors(n).end(); ++cl) {
        Cell* otherC = NULL;
        bool dtCell1st;

        if ((*dtCell) == (*cl)->first()) {
          otherC = (*cl)->second();
          dtCell1st = true;
        }
        else {
          otherC = (*cl)->first();
          dtCell1st = false;
        }
        assert(otherC);

    /*else { *//* We have different cells. */
    for (size_t c1 = 0; c1 < n_colours; ++c1) {
      for (size_t c2 = 0; c2 < n_colours; ++c2) {
	ColourPair *cp = M_MANAGER->cp(c1, c2); // We have pairs with same and different colours

	if (cp->needPairs()) {
	  double cutoff_sq = cp->cutoff() * cp->cutoff();

	  if (c1 < c2) {

            // createDistances for different cells cell1->free(c1) - cell2->free(c2)
            if (dtCell1st == true) {
              createDistances1stDirty
                (cp->freePairs(),
                cutoff_sq,
                1,
                (*cl)->first(),
                (*cl)->second(),
                (*cl)->first()->particles(c1),
                (*cl)->second()->particles(c2),
                (*cl)->actsOn().first,
                (*cl)->actsOn().second,
                (*cl)->mCellDist());
            }
            else {
              createDistances2ndDirty
                (cp->freePairs(),
                cutoff_sq,
                1,
                (*cl)->first(),
                (*cl)->second(),
                (*cl)->first()->particles(c1),
                (*cl)->second()->particles(c2),
                (*cl)->actsOn().first,
                (*cl)->actsOn().second,
                (*cl)->mCellDist());
            }

            // createDistances for different cells cell1->free(c1) - cell2->frozen(c2)
	    if ((*cl)->actsOn().first) {
              if (dtCell1st == true) {
                createDistances1stDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  1,
                  (*cl)->first(),
                  (*cl)->second(),
                  (*cl)->first()->particles(c1),
                  (*cl)->second()->frozenParticles(c2),
                  true,
                  false,
                  (*cl)->mCellDist());
              }
              else {
                createDistances2ndDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  1,
                  (*cl)->first(),
                  (*cl)->second(),
                  (*cl)->first()->particles(c1),
                  (*cl)->second()->frozenParticles(c2),
                  true,
                  false,
                  (*cl)->mCellDist());
              }
            } //END:  createDistances for different cells cell1->free(c1) - cell2->frozen(c2)

            // createDistances for different cells cell1->frozen(c1) - cell2->free(c2)
	    if ((*cl)->actsOn().second) {
              if (dtCell1st == true) {
                createDistances1stDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  1,
                  (*cl)->first(),
                  (*cl)->second(),
                  (*cl)->first()->frozenParticles(c1),
                  (*cl)->second()->particles(c2),
                  false,
                  true,
                  (*cl)->mCellDist());
              }
              else {
                createDistances2ndDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  1,
                  (*cl)->first(),
                  (*cl)->second(),
                  (*cl)->first()->frozenParticles(c1),
                  (*cl)->second()->particles(c2),
                  false,
                  true,
                  (*cl)->mCellDist());
              }
            } // END: createDistances for different cells cell1->frozen(c1) - cell2->free(c2)
	  } // END: if (c1 < c2)
          else {

              // createDistances for different cells cell2->free(c2) - cell1->free(c1)
//               list<Particle*> dirtyP = (*dtCell)->particles(c1);
//               list<Particle*> dirtyPFr2Cell = otherC->particles(c2);

              if (dtCell1st == true) {
                createDistances2ndDirty
                  (cp->freePairs(),
                  cutoff_sq,
                  -1,
                  (*cl)->second(),
                  (*cl)->first(),
                  (*cl)->second()->particles(c2),
                  (*cl)->first()->particles(c1),
                  (*cl)->actsOn().second,
                  (*cl)->actsOn().first,
                  (*cl)->mCellDist());
              }
              else {
                createDistances1stDirty
                  (cp->freePairs(),
                  cutoff_sq,
                  -1,
                  (*cl)->second(),
                  (*cl)->first(),
                  (*cl)->second()->particles(c2),
                  (*cl)->first()->particles(c1),
                  (*cl)->actsOn().second,
                  (*cl)->actsOn().first,
                  (*cl)->mCellDist());
              }




/*              list<Particle*>::iterator p1_end = dirtyPFr2Cell.end();
              list<Particle*>::iterator p2_end = dirtyP.end();

              for (list<Particle*>::iterator i = dirtyPFr2Cell.begin(); i != p1_end; ++i) {
                for (list<Particle*>::iterator j = dirtyP.begin(); j != p2_end; ++j) {
                  dist_t d;
                  d.abs_square = 0;
                  for (int _i = 0; _i < SPACE_DIMS; _i++) {
                    d.cartesian[_i] = (*cl)->mCellDist()[_i]
                      + (*i) -> r[_i] - otherC -> corner1[_i]
                      - (*j) -> r[_i] + (*dtCell) -> corner1[_i];
                    d.abs_square += d.cartesian[_i]*d.cartesian[_i];
                    }

                  if (d.abs_square < cutoff_sq) {
                        d.abs = sqrt(d.abs_square);
                        Pairdist* temp = &cp->freePairs(thread_no).newPair();
                        temp->set(d, *i, *j, (*cl)->actsOn().second, (*cl)->actsOn().first);

                        pairAddedFree(temp);
                  }
                }
              }*/ // END: createDistances for different cells cell2->free(c2) - cell1->free(c1)

            // createDistances for different cells cell2->frozen(c2) - cell1->free(c1)
	    if ((*cl)->actsOn().first) {
              list<Particle*> dirtyP2Cell = otherC->frozenParticles(c2);

              if (dtCell1st == true) {
                createDistances2ndDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  -1,
                  (*cl)->second(),
                  (*cl)->first(),
                  (*cl)->second()->frozenParticles(c2),
                  (*cl)->first()->particles(c1),
                  false,
                  true,
                  (*cl)->mCellDist());
              }
              else {
                createDistances1stDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  -1,
                  (*cl)->second(),
                  (*cl)->first(),
                  (*cl)->second()->frozenParticles(c2),
                  (*cl)->first()->particles(c1),
                  false,
                  true,
                  (*cl)->mCellDist());
              }
            }

            // createDistances for different cells cell2->free(c2) - cell1->frozen(c1)
	    if ((*cl)->actsOn().second) {
              if (dtCell1st == true) {
                createDistances2ndDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  -1,
                  (*cl)->second(),
                  (*cl)->first(),
                  (*cl)->second()->particles(c2),
                  (*cl)->first()->frozenParticles(c1),
                  true,
                  false,
                  (*cl)->mCellDist());
              }
              else {
                createDistances1stDirtyFrozen
                  (cp->frozenPairs(),
                  cutoff_sq,
                  -1,
                  (*cl)->second(),
                  (*cl)->first(),
                  (*cl)->second()->particles(c2),
                  (*cl)->first()->frozenParticles(c1),
                  true,
                  false,
                  (*cl)->mCellDist());
              }
            } // END: createDistances for different cells cell2->free(c2) - cell1->frozen(c1)
	  }
	}
      }
    }
  }
    } //Loop over 1st colour
      }
    }


time(&t4);

    // randomise pairs if wished
    if(M_PHASE->randomPairs())
    {
      FOR_EACH_COLOUR_PAIR
      (M_MANAGER,
//        MSG_DEBUG("ManagerCell::createDistances", "start, CP = " << cp->firstColour() << cp->secondColour());
       for (int t = 0; t < (int)global::n_threads; ++t)
       {
          // free pairs
          size_t listSize = cp->freePairs()[t].size();
          if(listSize)
          {
            cp->freePairsRandom(t).clear();
            for(size_t __i = 0; __i < listSize; ++__i)
              cp->freePairsRandom(t).newEntry().m_val = __i;
            for(size_t __i = listSize-1; __i > 0; --__i)
            {
              size_t slot = size_t(M_PHASE->rng().uniform()*(__i+1));
              size_t temp = cp->freePairsRandom(t)[__i].m_val;
              cp->freePairsRandom(t)[__i].m_val = cp->freePairsRandom(t)[slot].m_val;
              cp->freePairsRandom(t)[slot].m_val = temp;
            }
          }

          listSize = cp->frozenPairs()[t].size();
//           MSG_DEBUG("ManagerCell::createDistances", "listsize = " << listSize);
          if(listSize)
          {
            cp->frozenPairsRandom(t).clear();
            for(size_t __i = 0; __i < listSize; ++__i)
              cp->frozenPairsRandom(t).newEntry().m_val = __i;
            for(size_t __i = listSize-1; __i > 0; --__i)
            {
              size_t slot = size_t(M_PHASE->rng().uniform()*(__i+1));
              size_t temp = cp->frozenPairsRandom(t)[__i].m_val;
              cp->frozenPairsRandom(t)[__i].m_val = cp->frozenPairsRandom(t)[slot].m_val;
              cp->frozenPairsRandom(t)[slot].m_val = temp;
            }
          }
       }
      );
    }

    m_valid_dist = true;
    m_create_now = false;
  }
// time(&t7);
// time(&t1);
// MSG_INFO("VLYaoCreator", "creating new DISTS took " << difftime(t6, t7) << " seconds");
// MSG_INFO("VLYaoCreator", "RANDOMIZE took " << difftime(t1, t4) << " seconds");
// MSG_INFO("VLYaoCreator", "create Distances took " << difftime(t0, t1) << " seconds");
}

//  }
//}


void VLYaoCreator::updateAndSwitch(Particle* i, int& size /*i->tag.intByOffset(m_offset_free[i->c][c1])*//*, int& sizeP*/ /*i->tag.intByOffset(m_offset_free[c1][i->c])*/ , int size_loc /*m_offset_free[i->c][c1]*/, int size_locP /*m_offset_free[p->c][i->c]*/ , int capacityP /*m_capacity_free[p->c][i->c]*/, PairList* pl /*cp->freePairs(thread_no)*/) {

Cell *cell_one = M_MANAGER->findCell(i->r);

         for (int _i = 1; _i <= size; ++_i) {
	   Particle* p;
           int slot = i->tag.intByOffset(size_loc + _i*m_int_offset);

           Pairdist& pd = (*pl)[slot];

           int_point_t off;
           off.assign(0);
           int dir = 0;
           int n;
           dist_t d;
           d.abs_square = 0;
           int direction;
           bool same_cell = true;

// 		  Cell *cell_one = M_MANAGER->findCell(/*pair->*/pd.firstPart()->r);
           if (pd.firstPart() == i) {
             p = pd.secondPart();
             direction = 1;
           }
           else {
             p = pd.firstPart();
             direction = -1;
           }
//p is the second particle->the algorithm is right; p is the first particle->according to the settings above also should work.That is why instead of pd->secondPart() we can say "p".
           for (int j = 0; j < SPACE_DIMS; j++) {
             if (p/*d.secondPart()*/->r[j] < cell_one->corner1[j]) {
               if (p/*d.secondPart()*/->r[j] > (2*cell_one->corner1[j] - cell_one->corner2[j])) {
                 off[j] = -1;
               }
               else {
                 off[j] = 1;
               }
               same_cell = false;
             }
             else if (p/*d.secondPart()*/->r[j] >= cell_one->corner2[j]) {
               if (p/*d.secondPart()*/->r[j] <= (2*cell_one->corner2[j] - cell_one->corner1[j])) {
                 off[j] = 1;
               }
               else {
                 off[j] = -1;
               }
               same_cell = false;
             }
           }

           if (same_cell) {
             d.cartesian = pd.firstPart()->r - pd.secondPart()->r;
             for (int j = 0; j < SPACE_DIMS; j++)
               d.abs_square += d.cartesian[j]*d.cartesian[j];
           }
           else {
             OFFSET2NEIGHBOR(off, n);
 	     for (list<CellLink*>::iterator cl = cell_one->neighbors(n).begin(); cl != cell_one->neighbors(n).end(); ++cl) {
	     if(cell_one == (*cl)->first()) {
// 	       if((*cl)->second()->isInside(pd.secondPart()->r)/*was p here*/){
	       dir = 1;
	       for (int j = 0; j < SPACE_DIMS; j++) {
		 d.cartesian[j] = -dir*(*cl)->mCellDist()[j]*direction + pd.firstPart()->r[j] - pd.secondPart()->r[j] - (*cl)->first()->corner1[j]*direction + (*cl)->second()->corner1[j]*direction;
		 d.abs_square += d.cartesian[j]*d.cartesian[j];
		 cl = cell_one->neighbors(n).end();
		 --cl;
	       }
             }
	     else if (cell_one == (*cl)->second()) {
//             if ((*cl)->first()->isInside(pd.secondPart()->r)) {
               dir = -1;
               for (int j = 0; j < SPACE_DIMS; j++) {
                 d.cartesian[j] = -dir*(*cl)->mCellDist()[j]*direction + pd.firstPart()->r[j] - pd.secondPart()->r[j] - (*cl)->/*first*/second()->corner1[j]*direction + (*cl)->/*second*/first()->corner1[j]*direction;
                 d.abs_square += d.cartesian[j]*d.cartesian[j];
                 cl = cell_one->neighbors(n).end();
                 --cl;
               }
            }
              else
                throw gError("VLYaoCreator::updateAndSwitch", " The cell is not in this Celllink! Please contact the programmer.");
		  }
	        }
		  d.abs = sqrt(d.abs_square);
	  	  pd.set(d);

            int& sizeP = p->tag.intByOffset(size_locP);
            bool doSwitch = false;
            if (direction * d.cartesian.x > 0/*i->r.x > p->r.x*/) {
              doSwitch = true;

            }
            else if (d.cartesian.x == 0) {
              if (direction * d.cartesian.y > 0) {
                doSwitch = true;

              }
              else if (d.cartesian.y == 0) {
                if (direction * d.cartesian.z > 0) {
                  doSwitch = true;

                }
                else if (d.cartesian.z == 0)
                  throw gError("VLYaoCreator::createDistances Switching", " Positions of particle " + ObjToString(i->mySlot) + " and " +  ObjToString(p->mySlot) + " are identical. ");
              }
            }

            if (doSwitch) {
              ++sizeP/* += 1*/;
              if (sizeP < capacityP) {
                p->tag.intByOffset(size_locP + sizeP*m_int_offset) = slot;
              }
              else {
                throw gError("VLYaoCreator::createDistances", "There is not enough space to save all pairs of the current type. Please define a bigger slotSize!");
              }
              if (_i < size)
                i->tag.intByOffset(size_loc + _i*m_int_offset) = i->tag.intByOffset(size_loc + size*m_int_offset);
                --size;
                // We have to check the new element that comes to the same position
                --_i;
            }
         }
}


void VLYaoCreator::writeHeader()
{
    for (list<int>::iterator i = m_columns.begin(); i != m_columns.end(); i++)
      m_s << m_input_format->attrByIndex(*i).name << " ";
}



