/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#include "cbl_pair_particle_tensor.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()

const Callable_Register<CblPairParticleTensor> cbl_pair_particle_tensor("CblPairParticleTensor");


CblPairParticleTensor::CblPairParticleTensor(Simulation* parent)
  : CblPairParticleArbitrary(parent)
{
  m_function.setReturnType(Variant::TENSOR);
  m_1stparticleFactor.setReturnType(Variant::TENSOR);
  m_2ndparticleFactor.setReturnType(Variant::TENSOR);
  m_datatype = DataFormat::TENSOR;
  init();
}

CblPairParticleTensor::~CblPairParticleTensor()
{
}

void CblPairParticleTensor::init()
{
  m_properties.setClassName("CblPairParticleTensor");
  m_properties.setName("CblPairParticleTensor");

  m_properties.setDescription("Callable computing a user defined tensor property for particles, which has to be computed by pair summation. It may be defined by the attribute 'expression'. You should take care to define the attribute 'symmetry' properly. Additionally, notice that the particle expressions must be 3x3 tensors. If you only need scalars, use \"unitMat(scalar)\", where \"scalar\" is your scalar. The particle expressions are multiplied with the pair expression element by element.");

  m_1stPExpression = "unitMat(1)";
  m_2ndPExpression = "unitMat(1)";

  m_expression = "unitMat(1)";
}

void CblPairParticleTensor::setup()
{
  CblPairParticleArbitrary::setup();
}

