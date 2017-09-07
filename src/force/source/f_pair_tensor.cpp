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


#include "f_pair_tensor.h"

#include "threads.h"
#include "particle.h"
#include "simulation.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

const GenFTypeConcr<FPairTensor> f_pair_tensor("FPairTensor");

//---- Constructors/Destructor ----

FPairTensor::FPairTensor(Simulation *simulation): FPairArbitraryWF(simulation)
{
  init();
}

FPairTensor::~FPairTensor()
{
}

void FPairTensor::init()
{
  m_properties.setClassName("FPairTensor");

  m_properties.setDescription
      ("This is a completely general pair force FPV on a tensor t_i so that:\n"
      "dv_i = FPV*dt\n"
      "     = particleFactor_i*Sum_j(pairFactor_ij*weight_ij)*dt\n"
      "\nwhere particleFactor_i(j) is a tensor and a sum of quantities related to particle i(j),\n"     "pairFactor_ij is a tensor and includes all pair contributions of the pair ij,\n"
      "weight_ij represents -W'(rij)/rij of the used weighting function W.\n"
      );

  STRINGPC
      (tensor, m_tensor_name,
       "Name of the tensor field.");

//   FUNCTIONPAIRPC
//     (pairFactor, m_pairFactor,
//      "Function for pairFactor_ij. Type some nonsense to \n"
//      "obtain a complete list of possible variables and constants.\n"
//      "The expression may contain tensors and tensors,"
//      " but as a whole it must represent a tensor.");

  m_symmetry = 1;

  m_tensor_name = "tensor";

  m_pairFactorStr = "unitMat(1)";
  m_1stPExpression = "unitMat(1)";
  m_2ndPExpression = "unitMat(1)";
}


//---- Methods ----

#ifndef _OPENMP
void FPairTensor::computeForces(Pairdist* pair, int force_index)
#else
void FPairTensor::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  
  if (this->m_cutoff > pair->abs())
    {
      tensor_t temp;
      
      this->m_pairFactor(&temp, &(*pair));
      
      tensor_t fi;
      tensor_t fj;
      
      // compute the particle-expressions
      this->m_1stparticleFactor(&fi, &(*pair));
      this->m_2ndparticleFactor(&fj, &(*pair));
      
      fi *= temp;
      fj *= temp;
      
      fi *= this->m_wf->weight(pair, pair->secondPart()->r);
      fj *= this->m_wf->weight(pair, pair->firstPart()->r);
      
#ifndef _OPENMP
      if (pair->actsOnFirst())
	pair->firstPart()->tag.tensorByOffset(this->m_force_offset[force_index].first) += fi;
      if (pair->actsOnSecond())
	pair->secondPart()->tag.tensorByOffset(this->m_force_offset[force_index].second) += this->m_symmetry*fj;
#else
      if (pair->actsOnFirst()) {
	size_t _i = 0;
	for(size_t a = 0; a < SPACE_DIMS; ++a) {
	  for(size_t b = 0; b < SPACE_DIMS; ++b) {
	    (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi(a,b);
	    ++_i;
	  }
	}
      }
      if (pair->actsOnSecond()) {
	size_t _i = 0;
	for(size_t a = 0; a < SPACE_DIMS; ++a) {
	  for(size_t b = 0; b < SPACE_DIMS; ++b) {
	    (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += this->m_symmetry*fj(a,b);
	    ++_i;
	  }
	}
      }
#endif
      
    }
}


void FPairTensor::computeForces(Particle* part, int force_index)
{
  throw gError("FPairTensor::computeForces", "Fatal error: do not call FPairTensor::computeForces(Pairdist* pair, int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FPairTensor::computeForces(int force_index)
{
  throw gError("FPairTensor::computeForces", "Fatal error: do not call FPairTensor::computeForces(int force_index)!!! Needs a Particle argument. Please contact the programmer!");
}


void FPairTensor::setup()
{
  FPairArbitraryWF::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_pairFactor.setReturnType(Variant::TENSOR);
  m_1stparticleFactor.setReturnType(Variant::TENSOR);
  m_2ndparticleFactor.setReturnType(Variant::TENSOR);

  DataFormat::attribute_t firstAttr =
      Particle::s_tag_format[m_cp->firstColour()].attrByName(m_tensor_name);

  DataFormat::attribute_t secondAttr =
      Particle::s_tag_format[m_cp->secondColour()].attrByName(m_tensor_name);

  if(firstAttr.datatype != DataFormat::TENSOR)
    throw gError("FPairScalar::setup", "the symbol " + m_tensor_name +
        " is registerd as a non-tensor for species " +
            m_cp->manager()->species(m_cp->firstColour()));

  if(secondAttr.datatype != DataFormat::TENSOR)
    throw gError("FPairScalar::setup", "the symbol " + m_tensor_name +
        " is registerd as a non-tensor for species " +
            m_cp->manager()->species(m_cp->secondColour()));

  // FIXME: This is general for each FPairArbitrary except for FPairVels. That's why it's still here  => refine hierarchy or remove FPairVels
  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i].first =
        Particle::s_tag_format[m_cp->firstColour()].attrByName(string("force_"
        + m_tensor_name + "_" + ObjToString(i))).offset;
    m_force_offset[i].second =
        Particle::s_tag_format[m_cp->secondColour()].attrByName(string("force_"
        + m_tensor_name + "_" + ObjToString(i))).offset;
  }

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
//    MSG_DEBUG("FPairVels::setup", "Force colours same order as ColourPair's colours");
  }
  else {
    throw gError("FPairTensor::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FPairTensor::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (m_tensor_name == intr->dofIntegr()) {
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
        throw gError("FPairTensor::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
  if (col2 == intr->colour()) {
    if (m_tensor_name == intr->dofIntegr()) {
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
        throw gError("FPairTensor::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
}
#endif
