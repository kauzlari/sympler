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



#include "f_mdpd.h"

#include "threads.h"
#include "simulation.h"
#include "val_calculator_rho.h"
#include "pair_creator.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const GenFTypeConcr<FMDPD> f_generic("FMDPD");

//---- Constructors/Destructor ----

FMDPD::FMDPD(Simulation *simulation): FPairWF(simulation)
{
    init();
}


FMDPD::~FMDPD()
{
}


void FMDPD::init()
{
  m_properties.setClassName("FMDPD");

  m_properties.setDescription(
    "Implementation of MDPD."
  );

  FUNCTIONFIXEDPC
    (dpsi_dn, m_dpsi_dn,
     "dpsi/dn");

  DOUBLEPC
    (densityCutoff, m_density_cutoff, -1,
     "The density is not allowed to rise above this value.");

  BOOLPC
      (oneProp, m_oneProp,
       "Will the local density be computed over all ColourPairs (CP) or only over"
           " the CP corresponding to the chosen species?");
  STRINGPC
      (densitySymbol, m_rhoSymbol,
       "Name of the local density to be used.");


  m_density_cutoff = 10000.;
  m_oneProp = true;
  m_dpsi_dn.addVariable("n");
  m_rhoSymbol = "undefined";

  m_is_pair_force = true;
  m_is_particle_force = true;
}


//---- Methods ----

#ifndef _OPENMP
void FMDPD::computeForces(Pairdist* pair, int force_index)
#else
void FMDPD::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{
//   M_PAIRCREATOR->createDistances();
//
//   assert(m_cp);

  m_force_index = force_index;

//   FOR_EACH_PAIR__PARALLEL
//     (FMDPD,
//      m_cp,
     if (pair->abs() < this->m_cutoff) {
       point_t f;          /* Newtonian force */
       point_t fi;
       point_t fj;

       point_t rij = pair->cartesian();

       double weighti = this->m_wf->weight(pair, pair->firstPart()->r);
       double weightj = this->m_wf->weight(pair, pair->secondPart()->r);

       point_t ui = this->m_wf->localGradient(pair, pair->firstPart()->r);
       point_t uj = this->m_wf->localGradient(pair, pair->secondPart()->r);

       double ni = pair->firstPart()->tag.doubleByOffset(m_density_offset.first);
       if(ni > m_density_cutoff) ni = m_density_cutoff;
       double nj = pair->secondPart()->tag.doubleByOffset(m_density_offset.second);
       if(nj > m_density_cutoff) nj = m_density_cutoff;

       double dpsi_dni = this->m_dpsi_dn(ni);
       double dpsi_dnj = this->m_dpsi_dn(nj);

       f = (dpsi_dni*weightj+dpsi_dnj*weighti)*rij;
       fi = f + dpsi_dnj*ui;
       fj = -f + dpsi_dni*uj;

       fi = f;
       fj = -f;

       /* Set elements of the stress tensor */
//        for (int a = 0; a < SPACE_DIMS; ++a) {
//          for (int b = 0; b < SPACE_DIMS; ++b) {
// 	   pair->firstPart()->stress(a, b) += rij[a]*fi[b];
// 	   pair->secondPart()->stress(a, b) += rij[a]*fj[b];
//          }
//        }

       if (pair->actsOnFirst()) {
#ifndef _OPENMP
         pair->firstPart()->force[this->m_force_index] += fi;
#else
         for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
           (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi[_i];
         }
#endif
       }
       
       if (pair->actsOnSecond()) {
#ifndef _OPENMP
         pair->secondPart()->force[this->m_force_index] += fj;
#else
         for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
           (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += fj[_i];
         }
#endif
       }
     }
//     );

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
  assert(m_cp->firstColour() == m_cp->secondColour());


}


void FMDPD::computeForces(Particle* part, int force_index)
{
  /* We have to calculate the contribution of the particle itself! */
//   FOR_EACH_PARTICLE_C
//     (M_PHASE,
//      m_cp->firstColour(),

     point_t u = m_wf->localGradient(0, part->r);

     double n = part->tag.doubleByOffset(m_density_offset.first);
     if(n > m_density_cutoff) n = m_density_cutoff;

     double dpsi_dn = m_dpsi_dn(n);

     point_t f = dpsi_dn*u;

     part->force[m_force_index] += f;
//      );
}


void FMDPD::computeForces(int force_index)
{
  throw gError("FMDPD::computeForces", "Fatal error: do not call FMDPD::computeForces(int force_index)!!! Needs a Pairdist and a Particle argument. Please contact the programmer!");
}


void FMDPD::setup()
{
  FPairWF::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  if (m_dpsi_dn.isNull())
    throw gError(
      "FMDPD::setup",
      "Please define dpsi/dn.");

  // first we register the value in the CPs != m_cp, if m_oneProp = true
  if(m_oneProp)
  {
    for(size_t colour = 0; colour < M_MANAGER->nColours(); ++colour)
      if(!Particle::s_tag_format[colour].attrExists(m_rhoSymbol))
        throw gError("FMDPD::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + M_MANAGER->cp(colour, colour)->firstSpecies());


/*    FOR_EACH_COLOUR_PAIR
    (
      M_MANAGER,
      // if this is m_cp then do nothing (will be done afterwards)
      if(m_cp->firstColour() != cp->firstColour() || m_cp->secondColour() != cp->secondColour())
        cp->registerCalc(m_density_offset, new ValCalculatorRho(m_wf), true);
    );*/
  }
//   m_cp->registerCalc(m_density_offset, new ValCalculatorRho(m_wf), m_oneProp);
  else
  {
    if(!Particle::s_tag_format[m_cp->firstColour()].attrExists(m_rhoSymbol))
      throw gError("FMDPD::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + m_cp->firstSpecies());
    if(!Particle::s_tag_format[m_cp->secondColour()].attrExists(m_rhoSymbol))
      throw gError("FMDPD::setup", "Symbol '" + m_rhoSymbol + "' not found for species " + m_cp->secondSpecies());
  }
  m_density_offset.first = Particle::s_tag_format[m_cp->firstColour()].offsetByName(m_rhoSymbol);
  m_density_offset.second = Particle::s_tag_format[m_cp->secondColour()].offsetByName(m_rhoSymbol);


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
//    no change
  }
  else {
    throw gError("FMDPD::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FMDPD::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof = "vel_pos";
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
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
        throw gError("FMDPD::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");     
      }
    }
  }
  if (col2 == intr->colour()) {
    if (dof == intr->dofIntegr()) {
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
        throw gError("FMDPD::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");
      }
    }
  }
}
#endif


