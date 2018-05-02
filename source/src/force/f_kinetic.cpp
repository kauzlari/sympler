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



#include "f_kinetic.h"

#include "threads.h"
#include "simulation.h"
#include "val_calculator_rho.h"
#include "pair_creator.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const GenFTypeConcr<FKinetic> f_kinetic("FKinetic");

//---- Constructors/Destructor ----

FKinetic::FKinetic(Simulation *simulation): FPairWF(simulation)
{
  init();
}


FKinetic::~FKinetic()
{
}


void FKinetic::init()
{
  m_properties.setClassName("FKinetic");

  m_properties.setDescription(
    "Implementation of Boltzmann's kinetic equation."
  );

  STRINGPC
      (vvName, m_vv_name,
       "Name of the velocity correlation tensor.");

  STRINGPC
      (densitySymbol, m_rhoSymbol,
       "Name of the local density to be used.");

  DOUBLEPC
    (lambda, m_lambda, 0,
     "Relaxation constant for BGK approximation.");

  BOOLPC
      (oneProp, m_oneProp,
       "Will the local density be computed over all ColourPairs (CP) or only over"
           " the CP corresponding to the chosen species?");

  m_vv_name = "vv";
  m_lambda = 0.01;
  m_oneProp = true;
  m_rhoSymbol = "undefined";

  m_is_pair_force = true;
  m_is_particle_force = true;
}


//---- Methods ----

#ifndef _OPENMP
void FKinetic::computeForces(Pairdist* pair, int force_index)
#else
void FKinetic::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
//   M_PAIRCREATOR->createDistances();
//
//   assert(m_cp);

  m_force_index = force_index;

//   FOR_EACH_PAIR__PARALLEL
//     (FKinetic,
//      m_cp,

     if (pair->abs() < this->m_cutoff) {
       point_t fv;
       tensor_t fvv;

       fv.assign(0);
       fvv.assign(0);

       point_t rij = pair->cartesian();

       double weight = this->m_wf->weight(pair, fv);

       double ni = pair->firstPart()->tag.doubleByOffset(m_density_offset.first);
       double nj = pair->secondPart()->tag.doubleByOffset(m_density_offset.second);

       point_t vi = pair->firstPart()->v;
       point_t vj = pair->secondPart()->v;

       tensor_t vvi = pair->firstPart()->tag.tensorByOffset(m_vv_offset.first);
       tensor_t vvj = pair->secondPart()->tag.tensorByOffset(m_vv_offset.second);

       double kTmi = (vvi.trace() - vi.absSquare())/3;
       double kTmj = (vvj.trace() - vj.absSquare())/3;

       for (int a = 0; a < SPACE_DIMS; ++a) {
	 for (int b = 0; b < SPACE_DIMS; ++b) {
	   fv[a] +=
	     -weight*rij[b] * ((vvi(a, b) - vi[a]*vi[b])/nj - (vvj(a, b) - vj[a]*vj[b])/ni);

	   for (int c = 0; c < SPACE_DIMS; ++c) {
	     double vvvloci = vi[a]*vi[b]*vi[c];
	     double vvvlocj = vj[a]*vj[b]*vj[c];

	     if (b == c) {
	       vvvloci += kTmi*vi[a];
	       vvvlocj += kTmj*vj[a];
	     }
	     if (a == c) {
	       vvvloci += kTmi*vi[b];
	       vvvlocj += kTmj*vj[b];
	     }
	     if (a == b) {
	       vvvloci += kTmi*vi[c];
	       vvvlocj += kTmj*vj[c];
	     }

	     fvv(a, b) +=
	       -weight*rij[c] * ((vvvloci-vvi(a, b)*vi[c])/nj - (vvvlocj-vvj(a, b)*vj[c])/ni);
	   }
	 }
       }

       if (pair->actsOnFirst()) {
#ifndef _OPENMP
         pair->firstPart()->force[this->m_force_index] += fv;
	     pair->firstPart()->tag.tensorByOffset(this->m_vv_force_offset.first) += fvv;
#else
         for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
           (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fv[_i];
         }
         size_t __i = 0;
         for(size_t a = 0; a < SPACE_DIMS; ++a) {
           for(size_t b = 0; b < SPACE_DIMS; ++b) {
             (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + __i] += fvv(a, b);
             ++__i;
           }
         }
#endif
       }
       if (pair->actsOnSecond()) {
#ifndef _OPENMP
         pair->secondPart()->force[this->m_force_index] -= fv;
	     pair->secondPart()->tag.tensorByOffset(this->m_vv_force_offset.second) -= fvv;
#else
         for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
           (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] -= fv[_i];
         }
         size_t __i = 0;
         for(size_t a = 0; a < SPACE_DIMS; ++a) {
           for(size_t b = 0; b < SPACE_DIMS; ++b) {
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + __i] -= fvv(a, b);
             ++__i;
           }
         }
#endif
       }
     }
//      );

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  assert(m_cp->firstColour() == m_cp->secondColour());

}


void FKinetic::computeForces(Particle* part, int force_index)
{

     tensor_t fvv;

     point_t v = part->v;

     tensor_t vv = part->tag.tensorByOffset(m_vv_offset.first);

     double kTm = (vv.trace() - v.absSquare())/3;

     for (int a = 0; a < SPACE_DIMS; ++a) {
       for (int b = 0; b < SPACE_DIMS; ++b) {
	 fvv(a, b) = -m_lambda*vv(a, b);

	 if (a == b)
	   fvv(a, b) += m_lambda*v[a]*v[a]+kTm;
       }
     }

     part->tag.tensorByOffset(m_vv_force_offset.first) += fvv;
}


void FKinetic::computeForces(int force_index)
{
  throw gError("FKinetic::computeForces", "Fatal error: do not call FKinetic::computeForces(int force_index)!!! Needs a Pairdist and a Particle argument. Please contact the programmer!");
}


void FKinetic::setup()
{
  FPairWF::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  if(m_oneProp)
  {
    for(size_t colour = 0; colour < M_MANAGER->nColours(); ++colour)
      if(!Particle::s_tag_format[colour].attrExists(m_rhoSymbol))
        throw gError("FKinetic::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + M_MANAGER->cp(colour, colour)->firstSpecies());

  }
  else
  {
    if(!Particle::s_tag_format[m_cp->firstColour()].attrExists(m_rhoSymbol))
      throw gError("FKinetic::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + m_cp->firstSpecies());
    if(!Particle::s_tag_format[m_cp->secondColour()].attrExists(m_rhoSymbol))
      throw gError("FKinetic::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + m_cp->secondSpecies());

  }
  m_density_offset.first = Particle::s_tag_format[m_cp->firstColour()].offsetByName(m_rhoSymbol);
  m_density_offset.second = Particle::s_tag_format[m_cp->secondColour()].offsetByName(m_rhoSymbol);

  m_vv_offset.first
    = Particle::s_tag_format[m_cp->firstColour()].attrByName(m_vv_name).offset;
  m_vv_offset.second
    = Particle::s_tag_format[m_cp->secondColour()].attrByName(m_vv_name).offset;

  m_vv_force_offset.first
    = Particle::s_tag_format[m_cp->firstColour()].attrByName(STR_FORCE + STR_DELIMITER + m_vv_name).offset;
  m_vv_force_offset.second
    = Particle::s_tag_format[m_cp->secondColour()].attrByName(STR_FORCE + STR_DELIMITER + m_vv_name).offset;

// Setting the colours this force works on
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if ((col1 == cp->secondColour()) && (col2 == cp->firstColour())) {
    size_t dummy = col1;
    col1 = col2;
    col2 = dummy;
  }
  else if ((col1 == cp->firstColour()) && (col2 == cp->secondColour())) {
    // no changes
  }
  else {
    throw gError("FKinetic::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FKinetic::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof1 = "vel_pos";
  string dof2 = "m_vv_name";
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (dof1 == intr->dofIntegr()) {
      if (col1 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();
      }
      else if (col1 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();      
      }
      else {
        throw gError("FKinetic::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");     
      }
    }
    if (dof2 == intr->dofIntegr()) {
      if (col1 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();
      }
      else if (col1 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();      
      }
      else {
        throw gError("FKinetic::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");      
      }
    }
  }
  if (col2 == intr->colour()) {
    if (dof1 == intr->dofIntegr()) {
      if (col2 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();
      }
      else if (col2 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();      
      }
      else {
        throw gError("FKinetic::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");      
      }
    }
    if (dof2 == intr->dofIntegr()) {
      if (col2 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();
      }
      else if (col2 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();      
      }
      else {
        throw gError("FKinetic::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");      
      }
    }
  }
}
#endif


