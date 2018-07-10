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



#include "f_generic.h"

#include "threads.h"
#include "simulation.h"
#include "val_calculator_rho.h"
#include "pair_creator.h"

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER  M_SIMULATION->controller()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const GenFTypeConcr<FGeneric> f_generic("FGeneric");

//---- Constructors/Destructor ----

FGeneric::FGeneric(Simulation *simulation): FPairWF(simulation)
{
    init();
}


FGeneric::~FGeneric()
{
}


void FGeneric::init()
{
  m_properties.setClassName("FGeneric");

  m_properties.setDescription(
    "Implementation of MDPD together with DPDE. Additional correction terms are "
    "needed in that case. Right now, only a fixed noise amplitude is possible."
  );

  DOUBLEPC
    (densityCutoff, m_density_cutoff, -1,
     "The density is not allowed to rise above this value.");

  FUNCTIONFIXEDPC
    (de_dn, m_de_dn,
     "de/dn\nExample: For the van-der-Waals gas de/dn = -a*N, where N is the number of microscopic entities, the fluid particle consists of and a is the van-der-Waals parameter.");

  FUNCTIONFIXEDPC
    (Tds_dn, m_Tds_dn,
     "T*ds/dn = -Pressure\nExample: For the van-der-Waals gas T*ds/dn = -kB*T*N/(n*(1-n*b)), if we ignore long range contributions in the fluid particles' entropy. T is temperature, s is local entropy, n is local density, kB is the Boltzmann constant, N is the number of microscopic entities, the fluid particle consists of and b is the van-der-Waals parameter.");

  m_properties.addProperty
    ("surfaceEntropy", PropertyList::DOUBLE, &m_surface_entropy, NULL,
     "Entropy per surface density.");

  m_density_cutoff = 10000.;
  m_surface_entropy = 0;

  m_de_dn.addVariables("n", "e");
  m_Tds_dn.addVariables("n", "e");

  m_is_pair_force = true;
  m_is_particle_force = true;
}


//---- Methods ----

#ifndef _OPENMP
void FGeneric::computeForces(Pairdist* pair, int force_index)
#else
void FGeneric::computeForces(Pairdist* pair, int force_index, int thread_no)
#endif
{

     if (pair->abs() < this->m_cutoff) {
       point_t f;          /* Newtonian force */
       point_t fi;
       point_t fj;
       double fei;         /* Force on the internal energy */
       double fej;

       point_t rij = pair->cartesian();
       point_t vi = pair->firstPart()->v;
       point_t vj = pair->secondPart()->v;
       point_t vij = vi - vj;

       double weighti = this->m_wf->weight(pair, pair->firstPart()->r);
       double weightj = this->m_wf->weight(pair, pair->secondPart()->r);

       point_t ui = this->m_wf->localGradient(pair, pair->firstPart()->r);
       point_t uj = this->m_wf->localGradient(pair, pair->secondPart()->r);

       double Tdsdni;
       double Tdsdnj;
       double dAdni;
       double dAdnj;

       Tdsdni = this->m_pcee->Tds_dn(pair->firstPart());
       Tdsdnj = this->m_pcee->Tds_dn(pair->secondPart());
       dAdni = this->m_pcee->de_dn(pair->firstPart()) - Tdsdni;
       dAdnj = this->m_pcee->de_dn(pair->secondPart()) - Tdsdnj;

       f = (dAdni*weightj + dAdnj*weighti)*rij;
       fi = f + dAdnj*ui;
       fj = -f + dAdni*uj;

       fei = (rij*vij);
       fej = Tdsdnj * (weighti * fei + ui*vi);
       fei = Tdsdni * (weightj * fei + uj*vj);

       if (pair->actsOnFirst()) {
#ifndef _OPENMP
         pair->firstPart()->force[/*self->m_*/force_index] += fi;
	     pair->firstPart()->tag.doubleByOffset(this->m_eforce_offset[force_index].first) += fei;
#else
         for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
           (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + _i] += fi[_i];
         }
         (*pair->firstPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].first))[m_posInVec.first + SPACE_DIMS + 1] += fei;
#endif
       }
       if (pair->actsOnSecond()) {
#ifndef _OPENMP
         pair->secondPart()->force[/*self->m_*/force_index] += fj;
	     pair->secondPart()->tag.doubleByOffset(this->m_eforce_offset[force_index].second) += fej;
#else
         for(size_t _i = 0; _i < SPACE_DIMS; ++_i) {
           (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + _i] += fj[_i];
         }
         (*pair->secondPart()->tag.vectorDoubleByOffset(m_offsetToVec[thread_no].second))[m_posInVec.second + SPACE_DIMS + 1] += fej;
#endif
       }
     }

   ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));
   assert(m_cp->firstColour() == m_cp->secondColour());

}

// self-contribution
void FGeneric::computeForces(Particle* part, int force_index)
{

     double dist_from_wall;
     point_t u = m_wf->localGradient(NULL, part->r, &dist_from_wall);
     point_t sf;
     double T = 0;
     size_t wall_dir = 0;

     sf.assign(0);

     if (m_surface_entropy != 0) {
       wall_dir = dynamic_cast<WeightingFunctionWithWall*>(m_wf)->wallDir();

       sf = m_surface_entropy * m_wf->gradientSurfaceWeight(part->r);

       T = part->tag.doubleByOffset(m_temperature_offset.first);
     }

     double Tdsdn = m_pcee->Tds_dn(part);
     double dAdn = m_pcee->de_dn(part)-Tdsdn;

     point_t f = dAdn*u + T*sf;

     double fe = (Tdsdn*u + T*sf) * part->v;



     part->force[/*m_*/force_index] += f;
     part->tag.doubleByOffset(m_eforce_offset[force_index].first) += fe;

}


void FGeneric::computeForces(int force_index)
{
  throw gError("FGeneric::computeForces", "Fatal error: do not call FGeneric::computeForces(int force_index)!!! Needs a Pairdist and a Particle argument. Please contact the programmer!");
}


void FGeneric::setup()
{
  FPairWF::setup();

  ColourPair *m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second));

  if (m_de_dn.isNull())
    throw gError(
      "FGeneric::setup",
      "Please define de/dn.");
  if (m_Tds_dn.isNull())
    throw gError(
      "FGeneric::setup",
      "Please define T*ds/dn.");

  for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
  {
    m_eforce_offset[i].first =
        Particle::s_tag_format[m_cp->firstColour()].attrByName(string("force_internal_energy_" + ObjToString(i))).offset;
    m_eforce_offset[i].second =
        Particle::s_tag_format[m_cp->secondColour()].attrByName("force_internal_energy_" + ObjToString(i)).offset;
  }

  m_temperature_offset.first =
    Particle::s_tag_format[m_cp->firstColour()].attrByName("temperature").offset;
  m_temperature_offset.second =
    Particle::s_tag_format[m_cp->secondColour()].attrByName("temperature").offset;

  m_pcee = (ParticleCacheEnergyEntropy*) Particle::registerCache
    (new ParticleCacheEnergyEntropy
     (m_cp->firstColour(), m_density_cutoff, m_cp, m_wf, &m_de_dn, &m_Tds_dn));

// Setting the colours that this force works on
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
    throw gError("FGeneric::setup", "No ColourPair for these colours. Contact the programmer.");
  }
}


#ifdef _OPENMP
void FGeneric::setForceSlots(Integrator* intr, int thread_no) {
  size_t col1 = M_MANAGER->getColour(m_species.first);
  size_t col2 = M_MANAGER->getColour(m_species.second);
  string dof1 = "vel_pos";
  string dof2 = string("force_internal_energy_");
  ColourPair *cp = M_MANAGER->cp(col1, col2);

  if (col1 == intr->colour()) {
    if (dof1 == intr->dofIntegr()) {
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
        throw gError("FGeneric::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");        
      }
    }
    if (dof2 == intr->dofIntegr()) {
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
        throw gError("FGeneric::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");       
      }
    }
  }
  if (col2 == intr->colour()) {
    if (dof1 == intr->dofIntegr()) {
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
        throw gError("FGeneric::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");        
      }
    }
    if (dof2 == intr->dofIntegr()) {
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
        throw gError("FGeneric::setForceSlots", "No match for this Force's colours for the ColourPair. Contact the programmer.");        
      }
    }
  }
}
#endif



