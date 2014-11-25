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



#include "f_dpd.h"

#include "random.h"
#include "threads.h"
#include "simulation.h"
#include "val_calculator_rho.h"
#include "pair_creator.h"
#include "linked_list_creator.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const GenFTypeConcr<FDPD> f_dpd("FDPD");

//---- Constructors/Destructor ----

FDPD::FDPD(Simulation *simulation): FWithRng(simulation)
{
  init();
}


FDPD::~FDPD()
{
}


void FDPD::init()
{
  m_properties.setClassName("FDPD");

  m_properties.setDescription
    ("An implementation of the DPD forces."
     );

//   STRINGPC
//     (species1, m_species.first,
//      "First species, this force should act on.");
//
//   STRINGPC
//     (species2, m_species.second,
//      "Second species, this force should act on.");
//
//   STRINGPC
//     (weightingFunction, m_weighting_function,
//      "Defines the weighting function to be used.");

  DOUBLEPC
    (dissipation, m_dissipation, 0,
     "Defines the dissipation constant.");

  DOUBLEPC
    (kBToverM, m_temperature, 0,
     "k_BT/m = <v^2> to thermalize to.");

  m_weighting_function = "default";
  m_species.first = "fluid";
  m_species.second = "fluid";

  m_dissipation = 10;
  m_temperature = 1;
}


//---- Methods ----

#ifndef _OPENMP
void FDPD::computeForces(Pairdist* pair, int force_index)
#else
void FDPD::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
//   M_PAIRCREATOR->createDistances();
  //RandomNumberGenerator m_rng;
  //OLD:
//   assert(m_cp);

  m_force_index = force_index;

//   FOR_EACH_PAIR__PARALLEL
//     (FDPD,
//      m_cp,
     if (pair->abs() < this->m_cutoff) {
       /* We use interpolate because of the normalization */
       double weighti = this->m_wf->interpolate(pair, pair->firstPart()->r);
       double weightj = this->m_wf->interpolate(pair, pair->secondPart()->r);

       point_t vij = pair->firstPart()->v - pair->secondPart()->v;
       point_t eij = pair->cartesian()/pair->abs();

       double fd = -this->m_dissipation*(eij*vij);

       double fn = this->m_noise*m_rng.normal(1)*this->m_sqrt_r_dt;

       point_t fi = (fd*weightj+fn*sqrt(weightj))*eij;
       point_t fj = (fd*weighti+fn*sqrt(weighti))*eij;

       /* Set elements of the stress tensor */
//        for (int a = 0; a < SPACE_DIMS; ++a) {
//          for (int b = 0; b < SPACE_DIMS; ++b) {
// 	   pair->firstPart()->stress(a, b) += (*pair)[a]*fi[b];
// 	   pair->secondPart()->stress(a, b) += (*pair)[a]*fj[b];
//          }
//        }

#ifndef _OPENMP
       if (pair->actsOnFirst())
         pair->firstPart()->force[this->m_force_index] += fi;

       if (pair->actsOnSecond())
         pair->secondPart()->force[this->m_force_index] -= fj;
#else
       size_t c1 = M_MANAGER->getColour(m_species.first);
       size_t c2 = M_MANAGER->getColour(m_species.second);

       if (pair->actsOnFirst()) {
         for (int _i = 0; _i < SPACE_DIMS; ++_i) {
          if (c1 == pair->firstPart()->c)
            (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi[_i];
          else if (c2 == pair->firstPart()->c)
            (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += fi[_i];
         }
       }
       if (pair->actsOnSecond()) {
         for (int _i = 0; _i < SPACE_DIMS; ++_i) {
          if (c1 == pair->secondPart()->c)
            (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] -= fj[_i];
          else if (c2 == pair->secondPart()->c)
            (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] -= fj[_i];
         }
       }
#endif
     }
//     );
}


void FDPD::computeForces(Particle* part, int force_index)
{
  throw gError("FDPD::computeForces", "Fatal error: do not call FDPD::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FDPD::computeForces(int force_index)
{
  throw gError("FDPD::computeForces", "Fatal error: do not call FDPD::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FDPD::setup()
{
  GenF::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  m_wf = M_SIMULATION->findWeightingFunction(m_weighting_function);
  m_cutoff = m_wf->cutoff();

  /* Colour pair initialization */
  m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);


  m_cp->registerForce(this);
  m_cp->setNeedPairs(true);
  m_cp->setCutoff(m_cutoff);

  m_sqrt_r_dt = 1/sqrt(M_CONTROLLER->dt());

  m_noise = sqrt(2*m_temperature*m_dissipation);
}


#ifdef _OPENMP
void FDPD::setForceSlots(Integrator* intr, int thread_no) {
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
#endif


