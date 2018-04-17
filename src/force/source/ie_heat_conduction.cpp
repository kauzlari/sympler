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



#include "ie_heat_conduction.h"

#include "random.h"
#include "threads.h"
#include "particle.h"
#include "simulation.h"
#include "val_calculator_r_i.h"
#include "pair_creator.h"

#define INTERNAL_ENERGY "internal_energy"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

const GenFTypeConcr<IEHeatConduction> ie_heat_conduction("IEHeatConduction");

//---- Constructors/Destructor ----

IEHeatConduction::IEHeatConduction(Simulation *simulation)
  : FWithRng(simulation)
{
  init();
}


IEHeatConduction::~IEHeatConduction()
{
}


void IEHeatConduction::init()
{
  m_properties.setClassName("IEHeatConduction");

  m_properties.setDescription(
    "This object takes care of the heat conduction term "
    "of DPD with energy conservation."
  );

  FUNCTIONFIXEDPC
    (kappa, m_kappa,
     "Model for kappa(energy).");

  m_kappa.setExpression("1");
}


//---- Methods ----

#ifndef _OPENMP
void IEHeatConduction::computeForces(Pairdist* pair, int force_index)
#else
void IEHeatConduction::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
     if (pair->abs() < this->m_cutoff) {
       double weight =
         (pair->tag.doubleByOffset(this->m_compute_ri_offset) - this->m_rcinv)*pair->abs();

       double energy1 =
         pair->firstPart()->tag.doubleByOffset(this->m_energy_offset.first);
       double energy2 =
         pair->secondPart()->tag.doubleByOffset(this->m_energy_offset.second);

       double f = weight * weight * this->kappa(energy1, energy2)
         * (this->m_integrator.first->reciprocalTemperature(*pair->firstPart())
            - this->m_integrator.second->reciprocalTemperature(*pair->secondPart()))
         + weight * this->alpha(energy1, energy2) * m_rng.normal(1) * this->m_r_sqrt_dt;

       /* because of CONVENTION 2 in pairdist.h, actsOn*() will always return false for a
          frozen particle and there is consequently no danger, e.g., of modifying a force
          acting on a frozen particle, which does not exist in memory */

       if (pair->actsOnFirst()) {
#ifndef _OPENMP
         pair->firstPart()->tag.doubleByOffset(this->m_force_offset[force_index].first) += f;
#else
         size_t c1 = M_MANAGER->getColour(m_species.first);
         size_t c2 = M_MANAGER->getColour(m_species.second);

         if (c1 == pair->firstPart()->c)
           (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first] += f;
         else if (c2 == pair->firstPart()->c)
           (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second] += f;
#endif
       }
       if (pair->actsOnSecond()) {
#ifndef _OPENMP
         pair->secondPart()->tag.doubleByOffset(this->m_force_offset[force_index].second) -= f;
#else
         size_t c1 = M_MANAGER->getColour(m_species.first);
         size_t c2 = M_MANAGER->getColour(m_species.second);

         if (c1 == pair->secondPart()->c)
           (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first] -= f;
         else if (c2 == pair->secondPart()->c)
           (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second] -= f;
#endif
       }

     }
}

void IEHeatConduction::computeForces(Particle* part, int force_index)
{
  throw gError("IEHeatConduction::computeForces", "Fatal error: do not call IEHeatConduction::computeForces(Particle* pair, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void IEHeatConduction::computeForces(int force_index)
{
  throw gError("IEHeatConduction::computeForces", "Fatal error: do not call IEHeatConduction::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void IEHeatConduction::setup()
{
  FWithRng::setup();
  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_rcinv = 1/m_cutoff;
  m_r_sqrt_dt = 1/sqrt(M_CONTROLLER->dt());

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_force_offset[i].first =
        Particle::s_tag_format[m_cp->firstColour()].attrByName("force_internal_energy_"
          + ObjToString(i)).offset;
    m_force_offset[i].second =
        Particle::s_tag_format[m_cp->secondColour()].attrByName("force_internal_energy_"
          + ObjToString(i)).offset;
  }

  m_energy_offset.first =
    Particle::s_tag_format[m_cp->firstColour()].attrByName(INTERNAL_ENERGY).offset;
  m_energy_offset.second =
    Particle::s_tag_format[m_cp->secondColour()].attrByName(INTERNAL_ENERGY).offset;

  m_integrator.first = (IntegratorEnergy*)
    M_CONTROLLER->findIntegrator("IntegratorEnergy", m_cp->firstSpecies());
  m_integrator.second = (IntegratorEnergy*)
    M_CONTROLLER->findIntegrator("IntegratorEnergy", m_cp->secondSpecies());

  if (!m_integrator.first || !m_integrator.second)
    throw gError
	    ("IEHeatConduction",
	     "IEHeatConduction can only be used in connection with IntegratorEnergy.");

  m_cp->registerCalc(m_compute_ri_offset, new ValCalculatorRi, false);

  m_alpha.setExpression("sqrt(2*"+m_kappa.expression()+")");
}


#ifdef _OPENMP
void IEHeatConduction::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof = INTERNAL_ENERGY;

  if (col1 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
      intr->merge() = true;
      m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
      m_posInVec.first = intr->posInVec();
    }
  }
  if (col2 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
      intr->merge() = true;
      m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
      m_posInVec.second = intr->posInVec();
    }
  }
}
#endif

