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



#include "gen_f_solo.h"

#include "threads.h"
#include "simulation.h"
#include "val_calculator_r_i.h"
#include "pair_creator.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()


//---- Constructors/Destructor ----

GenFSolo::GenFSolo(Simulation *simulation): FWithRng(simulation), m_force_factors_old(true)
{
  init();
}


GenFSolo::~GenFSolo()
{
}


void GenFSolo::init()
{
  m_properties.setName("GenFSolo");
}


#ifndef _OPENMP
void GenFSolo::computeForces(Pairdist* pair, int force_index)
#else
void GenFSolo::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
  m_force_index = force_index;

     if (pair->abs() < this->m_cutoff) {
       point_t g = {{{ 0, 0, 0 }}};
       double weight = this->m_wf->interpolate(pair, g);

       point_t f = this->computeForceFactor(&(*pair)) * weight * pair->cartesian()/pair->abs();

       /* because of CONVENTION 2 in pairdist.h, actsOn*() will always return false for a
          frozen particle and there is consequently no danger, e.g., of modifying a force
          acting on a frozen particle, which does not exist in memory */
#ifndef _OPENMP
       if (pair->actsOnFirst())
         pair->firstPart()->force[this->m_force_index] += f;
       if (pair->actsOnSecond())
         pair->secondPart()->force[this->m_force_index] -= f;
#else

       // FIXME: next two lines unnecessary and slow !!!
       size_t c1 = M_MANAGER->getColour(m_species.first);
       size_t c2 = M_MANAGER->getColour(m_species.second);

       if (pair->actsOnFirst()) {
         for (int _i = 0; _i < SPACE_DIMS; ++_i) {
          if (c1 == pair->firstPart()->c)
            (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += f[_i];
          else if (c2 == pair->firstPart()->c)
            (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += f[_i];
         }
       }
       if (pair->actsOnSecond()) {
         for (int _i = 0; _i < SPACE_DIMS; ++_i) {
          if (c1 == pair->secondPart()->c)
            (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] -= f[_i];
          else if (c2 == pair->secondPart()->c)
            (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] -= f[_i];
         }
       }

#endif

     }
}


void GenFSolo::computeForces(Particle* part, int force_index)
{
  throw gError("GenFSolo::computeForces", "Fatal error: do not call GenFSolo::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void GenFSolo::computeForces(int force_index)
{
  throw gError("GenFSolo::computeForces", "Fatal error: do not call GenFSolo::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void GenFSolo::setup()
{
  FWithRng::setup();
  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

//   m_offset = m_cp->tagFormat().addAttribute
//     (string(FORCE_FACTOR_STR) + className() + "_" + m_cp->toString(), DataFormat::DOUBLE, false).offset;

  m_cp->registerCalc(m_compute_ri_offset, new ValCalculatorRi, false);

  m_rcinv = 1/m_cutoff;
}


#ifdef _OPENMP
void GenFSolo::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof = "vel_pos";

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


// void GenFSolo::mergeCopies(Particle* p, int thread_no, int force_index) {
//   size_t c1 = M_MANAGER->getColour(m_species.first);
//   size_t c2 = M_MANAGER->getColour(m_species.second);
//
//   for (int _i = 0; _i < SPACE_DIMS; ++_i) {
//     if (c1 == p->c) {
//       p->force[this->m_force_index][_i] += (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i];
//       (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] = 0;
//     }
//     else if (c2 == p->c) {
//       p->force[this->m_force_index][_i] += (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i];
//       (*p->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] = 0;
//     }
//   }
// }
#endif
