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



#include "f_pair_scalar.h"


#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "pair_creator.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

const GenFTypeConcr<FPairScalar> f_pair_scalar("FPairScalar");

//---- Constructors/Destructor ----

FPairScalar::FPairScalar(Simulation *simulation): FPairArbitraryWF(simulation)
{
  init();
}


FPairScalar::~FPairScalar()
{
}


void FPairScalar::init()
{
  m_properties.setClassName("FPairScalar");

  m_properties.setDescription(
    "This is a completely general pair force FPS on a scalar s_i so that:\n"
    "ds_i = FPS*dt\n"
    "     = particleFactor_i*Sum_j(pairFactor_ij*weight_ij)*dt\n"
      "\nwhere particleFactor_i(j) is a sum of quantities related to particle i(j),\n"
      "      pairFactor_ij includes all pair contributions of the pair ij,\n"
      "      weight_ij represents the derivative of the weighting function.\n"
    "\nNote that the expression has to give a scalar in the end,"
    "otherwise the compiler will report an error.");

  STRINGPC
    (scalar, m_scalar_name,
     "Name of the scalar field.");

  m_scalar_name = "scalar";

  m_1stPExpression = "1";
  m_2ndPExpression = "1";
  m_pairFactorStr = "1";
  m_symmetry = 1;

  m_is_pair_force = true;
  m_is_particle_force = false;
}


//---- Methods ----

#ifndef _OPENMP
void FPairScalar::computeForces(Pairdist* pair, int force_index)
#else
void FPairScalar::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
    if (this->m_cutoff > pair->abs())
    {
      double temp;

    double fi;
    double fj;

    // compute the particle-expressions
    this->m_1stparticleFactor(&fi, &(*pair));
    this->m_2ndparticleFactor(&fj, &(*pair));


    this->m_pairFactor(&temp, &(*pair));

    fi *= temp*this->m_wf->weight(pair, pair->secondPart()->r);
    fj *= temp*this->m_wf->weight(pair, pair->firstPart()->r);

#ifndef _OPENMP
      if (pair->actsOnFirst())
        pair->firstPart()->tag.doubleByOffset(this->m_force_offset[force_index].first) += fi;

      if (pair->actsOnSecond())
        pair->secondPart()->tag.doubleByOffset(this->m_force_offset[force_index].second) += this->m_symmetry*fj;

#else
       if (pair->actsOnFirst()) {
           (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first] += fi;
       }
       if (pair->actsOnSecond()) {
           (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second] += this->m_symmetry*fj;
       }
#endif

     }
}


void FPairScalar::computeForces(Particle* part, int force_index)
{
  throw gError("FPairScalar::computeForces", "Fatal error: do not call FPairScalar::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairScalar::computeForces(int force_index)
{
  throw gError("FPairScalar::computeForces", "Fatal error: do not call FPairScalar::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FPairScalar::setup()
{
  FPairArbitraryWF::setup();

  ColourPair *cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_1stparticleFactor.setReturnType(Variant::SCALAR);
  m_2ndparticleFactor.setReturnType(Variant::SCALAR);
  m_pairFactor.setReturnType(Variant::SCALAR);

  DataFormat::attribute_t firstAttr =
    Particle::s_tag_format[cp->firstColour()].attrByName(m_scalar_name);

  DataFormat::attribute_t secondAttr =
    Particle::s_tag_format[cp->secondColour()].attrByName(m_scalar_name);

  if(firstAttr.datatype != DataFormat::DOUBLE)
    throw gError("FPairScalar::setup", "the symbol " + m_scalar_name +
    " is registerd as a non-scalar for species " +
    cp->manager()->species(cp->firstColour()));

  if(secondAttr.datatype != DataFormat::DOUBLE)
    throw gError("FPairScalar::setup", "the symbol " + m_scalar_name +
    " is registerd as a non-scalar for species " +
    cp->manager()->species(cp->secondColour()));

  // FIXME: This is general for each FPairArbitrary except for FPairVels. That's why it's still here  => refine hierarchy or remove FPairVels
  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i].first =
      Particle::s_tag_format[cp->firstColour()].attrByName(string("force_"
        + m_scalar_name + "_" + ObjToString(i))).offset;
    m_force_offset[i].second =
      Particle::s_tag_format[cp->secondColour()].attrByName(string("force_"
        + m_scalar_name + "_" + ObjToString(i))).offset;
  }
  
}


#ifdef _OPENMP
void FPairScalar::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
// MSG_DEBUG("FPairScalar::setForceSlots", " col1 = col integr c = " << col1 << " intr = " << intr->name());
    if (m_scalar_name == intr->dofIntegr()) {
// MSG_DEBUG("FPairScalar::setForceSlots", " scalar name = dof integr by col1. scalar = " << m_scalar_name);
      if (col1 == cp->firstColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
        m_posInVec.first = intr->posInVec();
// MSG_DEBUG("FPairScalar::setForceSlots", "col1. m_posInVec = " << m_posInVec.first);
      }
      else if (col1 == cp->secondColour()) {
        intr->merge() = true;
        m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
        m_posInVec.second = intr->posInVec();        
      }
      else {
        throw gError("FPairScalar::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");        
      }
    }
  }
  if (col2 == intr->colour()) {
    if (m_scalar_name == intr->dofIntegr()) {
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
        throw gError("FPairScalar::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");        
      }
    }
  }
}
#endif
