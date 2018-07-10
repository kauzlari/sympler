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


#include "cbl_pair_particle_vector.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()

const Callable_Register<CblPairParticleVector> cbl_pair_particle_vector("CblPairParticleVector");


CblPairParticleVector::CblPairParticleVector(Simulation* parent)
  : CblPairParticleArbitrary(parent)
{
  m_function.setReturnType(Variant::VECTOR);
  m_1stparticleFactor.setReturnType(Variant::VECTOR);
  m_2ndparticleFactor.setReturnType(Variant::VECTOR);
  m_datatype = DataFormat::POINT;
  init();
}

CblPairParticleVector::~CblPairParticleVector()
{
}

void CblPairParticleVector::init()
{
  m_properties.setClassName("CblPairParticleVector");
  m_properties.setName("CblPairParticleVector");

  m_properties.setDescription("Callable computing a user defined vector property for particles, which has to be computed by pair summation. It may be defined by the attribute 'expression'. You should take care to define the attribute 'symmetry' properly. Additionally, notice that the particle expressions must be vectors. If you only need scalars, use \"idVec(scalar)\", where \"scalar\" is your scalar. The particle expressions are multiplied with the pair expression element by element.");

  m_1stPExpression = "idVec(1)";
  m_2ndPExpression = "idVec(1)";

  m_expression = "idVec(1)";
}

void CblPairParticleVector::setup()
{
  CblPairParticleArbitrary::setup();
}

