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



#include "f_pair_vector.h"

#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "pair_creator.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

const GenFTypeConcr<FPairVector> f_pair_vector("FPairVector");

//---- Constructors/Destructor ----

FPairVector::FPairVector(Simulation *simulation): FPairArbitraryWF(simulation)
{
  init();
}

FPairVector::~FPairVector()
{
}

void FPairVector::init()
{
  m_properties.setClassName("FPairVector");

  m_properties.setDescription
    ("This is a completely general pair force FPV on a vector v_i so that:\n"
     "dv_i = FPV*dt\n"
     "     = particleFactor_i*Sum_j(pairFactor_ij*weight_ij)*dt\n"
      "\nwhere particleFactor_i(j) is a vector and a sum of quantities related to particle i(j),\n"     "pairFactor_ij is a vector and includes all pair contributions of the pair ij,\n"
     "weight_ij represents the derivative of the used weighting function.\n"
     );

  STRINGPC
    (vector, m_vector_name,
     "Name of the vector field.");

//   FUNCTIONPAIRPC
//     (pairFactor, m_pairFactor,
//      "Function for pairFactor_ij. Type some nonsense to \n"
//      "obtain a complete list of possible variables and constants.\n"
//      "The expression may contain vectors and tensors,"
//      " but as a whole it must represent a vector.");

  m_symmetry = 1;

  m_vector_name = "vector";

  m_pairFactorStr = "idVec(0)";
  m_1stPExpression = "idVec(0)";
  m_2ndPExpression = "idVec(0)";

  m_is_pair_force = true;
  m_is_particle_force = false;
}


//---- Methods ----

#ifndef _OPENMP
void FPairVector::computeForces(Pairdist* pair, int force_index)
#else
void FPairVector::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
       if (this->m_cutoff > pair->abs())
       {
         point_t temp;

    this->m_pairFactor(&temp, &(*pair));

    point_t fi;
    point_t fj;

    // compute the particle-expressions
    this->m_1stparticleFactor(&fi, &(*pair));
    this->m_2ndparticleFactor(&fj, &(*pair));

    // loop necessary because operator* of math_vector_t does scalar product
    for(size_t i = 0; i < SPACE_DIMS; ++i)
    {
      fi[i] *= temp[i];
      fj[i] *= temp[i];
    }

    fi *= this->m_wf->weight(pair, pair->secondPart()->r);
    fj *= this->m_wf->weight(pair, pair->firstPart()->r);

#ifndef _OPENMP
    if (pair->actsOnFirst())
      pair->firstPart()->tag.pointByOffset(this->m_force_offset[force_index].first) += fi;
    if (pair->actsOnSecond())
      pair->secondPart()->tag.pointByOffset(this->m_force_offset[force_index].second) += this->m_symmetry*fj;
#else
    if (pair->actsOnFirst()) {
      for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
        (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi[_i];
      }
    }
    if (pair->actsOnSecond()) {
      for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
        (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += this->m_symmetry*fj[_i];
      }
    }
#endif

       }
}


void FPairVector::computeForces(Particle* part, int force_index)
{
  throw gError("FPairVector::computeForces", "Fatal error: do not call FPairVector::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairVector::computeForces(int force_index)
{
  throw gError("FPairVector::computeForces", "Fatal error: do not call FPairVector::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairVector::setup()
{
  FPairArbitraryWF::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_pairFactor.setReturnType(Variant::VECTOR);
//   m_pairFactor.setColourPair(m_cp);
  m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.setReturnType(Variant::VECTOR);


/*  if (m_pairFactor.isNull()) {
    throw gError("FPairVector::setup", "Please specify a function for 'pairFactor'.");
  }*/

  DataFormat::attribute_t firstAttr =
    Particle::s_tag_format[m_cp->firstColour()].attrByName(m_vector_name);

  DataFormat::attribute_t secondAttr =
    Particle::s_tag_format[m_cp->secondColour()].attrByName(m_vector_name);

  if(firstAttr.datatype != DataFormat::POINT)
    throw gError("FPairScalar::setup", "the symbol " + m_vector_name +
    " is registerd as a non-vector for species " +
    m_cp->manager()->species(m_cp->firstColour()));

  if(secondAttr.datatype != DataFormat::POINT)
    throw gError("FPairScalar::setup", "the symbol " + m_vector_name +
    " is registerd as a non-vector for species " +
    m_cp->manager()->species(m_cp->secondColour()));

  // FIXME: This is general for each FPairArbitrary except for FPairVels. That's why it's still here  => refine hierarchy or remove FPairVels
  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i].first =
      Particle::s_tag_format[m_cp->firstColour()].attrByName(string("force_"
        + m_vector_name + "_" + ObjToString(i))).offset;
    m_force_offset[i].second =
      Particle::s_tag_format[m_cp->secondColour()].attrByName(string("force_"
        + m_vector_name + "_" + ObjToString(i))).offset;
  }

// Setting the colours this force works on
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if ((col1 == cp->secondColour()) && (col2 == cp->firstColour())) {
    size_t dummy = col1;
    col1 = col2;
    col1 = dummy;
  }
  else if ((col1 == cp->firstColour()) && (col2 == cp->secondColour())) {
//    MSG_DEBUG("FPairVels::setup", "Force colours same order as ColourPair's colours");
  }
  else {
    throw gError("FPairVector::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FPairVector::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (m_vector_name == intr->dofIntegr()) {
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
        throw gError("FPairVector::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
  if (col2 == intr->colour()) {
    if (m_vector_name == intr->dofIntegr()) {
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
        throw gError("FPairVector::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
}



// void FPairVector::mergeCopies(Particle* p, int thread_no, int force_index) {
//   size_t c1 = M_MANAGER->getColour(m_species.first);
//   size_t c2 = M_MANAGER->getColour(m_species.second);
// //FIXME: .first and .second at the left side
//   for (int _i = 0; _i < SPACE_DIMS; ++_i) {
//     if (c1 == p->c) {
//       p->tag.pointByOffset(this->m_force_offset[force_index].first)[_i] += (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i];
//       (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] = 0;
//     }
//     else if (c2 == p->c) {
//       p->tag.pointByOffset(this->m_force_offset[force_index].second)[_i] += (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i];
//       (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] = 0;
//     }
//   }
// }
#endif

