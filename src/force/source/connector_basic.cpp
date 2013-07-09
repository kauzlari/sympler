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



#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "function_pair.h"

#include "connector_basic.h"



#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

const GenFTypeConcr<ConnectBasic> connector_basic("ConnectBasic");

//---- Constructors/Destructor ----

ConnectBasic::ConnectBasic()
{
  throw gError("ConnectorBasic::ConnectorBasic(default)", "Should not be called. Contact the programmer.");
}

ConnectBasic::ConnectBasic(Simulation *simulation): GenConnector(simulation)
{
  init();
}


ConnectBasic::~ConnectBasic()
{
}


void ConnectBasic::init()
{
  m_properties.setClassName("ConnectBasic");

  m_properties.setDescription(
    "This is a general basic connection force.");

  FUNCTIONPAIRPC (pairFactor, m_pairFactor,
     "Function for pairFactor_ij. Type some nonsense to \n"
     "obtain a complete list of possible variables and constants.\n"
     "The expression may contain vectors and tensors,"
     " but as a whole it must represent a vector.");

m_is_pair_force = false;
m_is_particle_force = false;
}


//---- Methods ----

/*FIXME: when parallelised, this function gets a Pairdist as argument and doesn't loop over a PairList by itself. For now (before first working parallel revision) we have to keep it as it is.*/
#ifndef _OPENMP
void ConnectBasic::computeForces(Pairdist* pair, int force_index)
#else
void ConnectBasic::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  throw gError("ConnectBasic::computeForces", "Fatal error: do not call ConnectBasic::computeForces(Pairdist* pair, int force_index)!!! Please contact the programmer!");
}

// OLD style force before parallelisation; still active
void ConnectBasic::computeForces(int force_index)
{
  // FIXME: not so nice to do that here, but currently the only possibility. See also GenConnector::setupAfterParticleCreation(). Is the problem solved if this function will be called for each pair (when parallelised)?
  if(!m_connectedList)
    m_connectedList = m_cp->connectedList(m_force_name);
  
  for (Pairdist *pair = m_connectedList->first(); pair != NULL; pair = pair->next)
  {
	point_t temp;
	m_pairFactor(&temp, &(*pair));
	point_t fi = temp;
	point_t fj = temp;
	pair->firstPart()->force[force_index] += fi;
	pair->secondPart()->force[force_index] -= fj;
  }
}


void ConnectBasic::computeForces(Particle* part, int force_index)
{
  throw gError("ConnectBasic::computeForces", "Fatal error: do not call ConnectBasic::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void ConnectBasic::setup()
{
  GenConnector::setup();

  m_pairFactor.setReturnType(Variant::VECTOR);
  m_pairFactor.setColourPair(m_cp);

  if (m_pairFactor.isNull()) {
    throw gError("ConnectorBasic::setup", "Please specify a function for 'pairFactor'.");
  }
}


#ifdef _OPENMP
void ConnectBasic::setForceSlots(Integrator* intr, int thread_no) {
//   size_t col1 = M_MANAGER->getColour(m_species.first);
//   size_t col2 = M_MANAGER->getColour(m_species.second);
//   string dof = "vel_pos";
//
//   if (col1 == intr->colour()) {
//     if (dof == intr->dofIntegr()) {
//       m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
//       m_posInVec = intr->posInVec();
//     }
//   }
//   if (col2 == intr->colour()) {
//     if (dof == intr->dofIntegr()) {
//       m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
//       m_posInVec = intr->posInVec();
//     }
//   }
}
#endif
