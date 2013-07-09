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




#include "f_dpde.h"

#include "random.h"
#include "threads.h"
#include "simulation.h"
#include "pair_creator.h"
// #include "val_calculator_rho.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const GenFTypeConcr<FDPDE> f_dpde("FDPDE");

//---- Constructors/Destructor ----

FDPDE::FDPDE(Simulation *simulation): FWithRng(simulation)
{
  init();
}


FDPDE::~FDPDE()
{
}


void FDPDE::init()
{
  m_properties.setClassName("FDPDE");

  m_properties.setDescription
    ("An implementation of the DPD with energy conservation forces."
     );

  DOUBLEPC
    (dissipation, m_dissipation, 0,
     "Defines the dissipation constant.");

  DOUBLEPC
    (noise, m_noise, 0,
     "Defines the noise amplitude.");

  m_dissipation = -1;
  m_noise = -1;
}


//---- Methods ----

#ifndef _OPENMP
void FDPDE::computeForces(Pairdist* pair, int force_index)
#else
void FDPDE::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
//   M_PAIRCREATOR->createDistances();
  //RandomNumberGenerator m_rng;

  //OLD:
//   assert(m_cp);

  m_force_index = force_index;

  if (m_noise < 0) {
    MSG_DEBUG("FDPDE::computeForces", "Fixed dissipation.");

//     FOR_EACH_PAIR__PARALLEL
//       (FDPDE,
//        m_cp,
       if (pair->abs() < this->m_cutoff) {
	 point_t g;
	 g.assign(0);


	 double weight = this->m_wf->interpolate(pair, g);

	 point_t vij = pair->firstPart()->v - pair->secondPart()->v;
	 point_t eij = pair->cartesian()/pair->abs();

	 double ev = eij*vij;

	 double noise_sq =
	   4*this->m_dissipation/
	   (this->m_ie.first->reciprocalTemperature(*pair->firstPart())
	    +this->m_ie.second->reciprocalTemperature(*pair->secondPart())
	    );

	 double fd = -this->m_dissipation*weight*ev;

	 double fn = m_rng.normal(1)*sqrt(noise_sq*weight)*this->m_sqrt_r_dt;

	 double fe = 0.5*
	   (-fd*ev
	    -noise_sq*weight
	    -fn*ev);

	 /* Fixme!!! Calculate stress tensor! */
	 /* Fixme!!! Take care of walls */

	 if (pair->actsOnFirst()) {
#ifndef _OPENMP
	   pair->firstPart()->force[this->m_force_index] += (fd+fn)*eij;
	   pair->firstPart()->tag.doubleByOffset(this->m_eforce_offset.first) += fe;
#else
           size_t c1 = M_MANAGER->getColour(m_species.first);
           size_t c2 = M_MANAGER->getColour(m_species.second);

             if (c1 == pair->firstPart()->c) {
               for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
                 (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += (fd+fn)*eij[_i];
               }
               (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + SPACE_DIMS + 1] += fe;
             }
             else if (c2 == pair->firstPart()->c) {
               for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
                 (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += (fd+fn)*eij[_i];
               }
               (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + SPACE_DIMS + 1] += fe;
             }
#endif
	 }

	 if (pair->actsOnSecond()) {
#ifndef _OPENMP
	   pair->secondPart()->force[this->m_force_index] -= (fd+fn)*eij;
	   pair->secondPart()->tag.doubleByOffset(this->m_eforce_offset.second) += fe;
#else
           size_t c1 = M_MANAGER->getColour(m_species.first);
           size_t c2 = M_MANAGER->getColour(m_species.second);

           if (c1 == pair->secondPart()->c) {
             for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] -= (fd+fn)*eij[_i];
           }
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + SPACE_DIMS + 1] += fe;
           }
           else if (c2 == pair->secondPart()->c) {
             for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] -= (fd+fn)*eij[_i];
           }
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + SPACE_DIMS + 1] += fe;
           }
#endif
	 }
       }
//        );
  } else {
    MSG_DEBUG("FDPDE::computeForces", "Fixed noise.");

//     FOR_EACH_PAIR__PARALLEL
//       (FDPDE,
//        m_cp,
       if (pair->abs() < this->m_cutoff) {
	 point_t g;
	 g.assign(0);

	 double weight = this->m_wf->interpolate(pair, g);

	 point_t vij = pair->firstPart()->v - pair->secondPart()->v;
	 point_t eij = pair->cartesian()/pair->abs();

	 double ev = eij*vij;

	 double dissipation =
	   this->m_noise*this->m_noise*
	   (this->m_ie.first->reciprocalTemperature(*pair->firstPart())
	    +this->m_ie.second->reciprocalTemperature(*pair->secondPart())
	    )/4;

	 double fd = -dissipation*weight*ev;

	 double fn = m_noise*m_rng.normal(1)*sqrt(weight)*this->m_sqrt_r_dt;

	 double fe = 0.5*
	   (-fd*ev
	    -m_noise*m_noise*weight
	    -fn*ev);

	 point_t f = (fd+fn)*eij;

	 /* Set elements of the stress tensor */
// 	 for (int a = 0; a < SPACE_DIMS; ++a) {
// 	   for (int b = 0; b < SPACE_DIMS; ++b) {
// 	     double h = (*pair)[a]*f[b];
// 	     pair->firstPart()->stress(a, b) += h;
// 	     pair->secondPart()->stress(a, b) += h;
// 	   }
// 	 }

	 if (pair->actsOnFirst()) {
#ifndef _OPENMP
	   pair->firstPart()->force[this->m_force_index] += f;
	   pair->firstPart()->tag.doubleByOffset(this->m_eforce_offset.first) += fe;
#else
           size_t c1 = M_MANAGER->getColour(m_species.first);
           size_t c2 = M_MANAGER->getColour(m_species.second);

             if (c1 == pair->firstPart()->c) {
               for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
                 (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += f[_i];
               }
               (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + SPACE_DIMS + 1] += fe;
             }
             else if (c2 == pair->firstPart()->c) {
               for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
                 (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += f[_i];
               }
               (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + SPACE_DIMS + 1] += fe;
             }
#endif
	 }

	 if (pair->actsOnSecond()) {
#ifndef _OPENMP
	   pair->secondPart()->force[this->m_force_index] -= f;
	   pair->firstPart()->tag.doubleByOffset(this->m_eforce_offset.first) += fe;
#else
           size_t c1 = M_MANAGER->getColour(m_species.first);
           size_t c2 = M_MANAGER->getColour(m_species.second);

           if (c1 == pair->secondPart()->c) {
             for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] -= f[_i];
           }
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + SPACE_DIMS + 1] += fe;
           }
           else if (c2 == pair->secondPart()->c) {
             for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] -= f[_i];
           }
             (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + SPACE_DIMS + 1] += fe;
           }
#endif
	 }
       }
//        );
  }
}

void FDPDE::computeForces(Particle* part, int force_index)
{
  throw gError("FDPDE::computeForces", "Fatal error: do not call FDPDE::computeForces(Particle* part, int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FDPDE::computeForces(int force_index)
{
  throw gError("FDPDE::computeForces", "Fatal error: do not call FDPDE::computeForces(int force_index)!!! Needs a Pairdist argument. Please contact the programmer!");
}


void FDPDE::setup()
{
  FWithRng::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  /* Get integrators */
  m_ie.first =
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species.first);
  m_ie.second =
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species.second);

  /* Get degrees of freedom */
  m_eforce_offset.first =
    Particle::s_tag_format[m_cp->firstColour()].attrByName("force_internal_energy").offset;
  m_eforce_offset.second =
    Particle::s_tag_format[m_cp->secondColour()].attrByName("force_internal_energy").offset;

  m_sqrt_r_dt = 1/sqrt(M_CONTROLLER->dt());

  if (m_dissipation < 0 && m_noise < 0)
    throw gError
      ("FDPDE::setup",
       "Please specify a dissipation constant or a noise amplitude.");

  if (m_dissipation > 0 && m_noise > 0)
    throw gError
      ("FDPDE::setup",
       "Please specify either a dissipation constant or a noise amplitude, not both.");
}


#ifdef _OPENMP
void FDPDE::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof1 = "vel_pos";
  string dof2 = "force_internal_energy";

  if (col1 == intr->colour()) {
    if (dof1 == intr->dofIntegr()) {
      intr->merge() = true;
      m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
      m_posInVec.first = intr->posInVec();
    }
    if (dof2 == intr->dofIntegr()) {
      intr->merge() = true;
      m_offsetToVec[thread_no].first = intr->offsetToVec()[thread_no];
      m_posInVec.first = intr->posInVec();
    }
  }
  if (col2 == intr->colour()) {
    if (dof1 == intr->dofIntegr()) {
      intr->merge() = true;
      m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
      m_posInVec.second = intr->posInVec();
    }
    if (dof2 == intr->dofIntegr()) {
      intr->merge() = true;
      m_offsetToVec[thread_no].second = intr->offsetToVec()[thread_no];
      m_posInVec.second = intr->posInVec();
    }
  }
}
#endif



