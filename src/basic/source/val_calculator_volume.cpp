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


#include "val_calculator_volume.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"

const SymbolRegister<ValCalculatorVolume> val_calc_volume("ValCalculatorVolume");

#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PAIRCREATOR M_PHASE->pairCreator()

// the order of the arguments was changed by purpose to detect the calls that must be modified
ValCalculatorVolume::ValCalculatorVolume(string densitySymbol/*extra_string*/, string volumeSymbol, WeightingFunction *wf)
  : NonBondedPairParticleCalculator(volumeSymbol), m_cutoff(wf->cutoff()), m_wf(wf), m_rhoSymbol(densitySymbol)/*, m_symbolName*//*m_extra_string*//*(densitySymbol*//*extra_string*//*)*/
{
  m_stage = 1;
  MSG_DEBUG("ValCalculatorVolume::ValCalculatorVolume", "CONSTRUCTOR");
}

ValCalculatorVolume::ValCalculatorVolume(/*Node*/Simulation* parent)
  : NonBondedPairParticleCalculator(parent)
{
  m_stage = 1;
  m_datatype = DataFormat::DOUBLE;
  init();
}

void ValCalculatorVolume::init()
{
//   m_properties.setClassName("ValCalculatorVolume");

  m_properties.setClassName("ValCalculatorPart");
  m_properties.setName("ValCalculatorVolume");

  m_properties.setDescription
    ("Calculates the local volume V_i according to the formula\n"
     "\n"
     "V_i=sum_j(W_ij/n_j)\n"
     "\n"
     "where W_ij is the value for particle pair (i,j) of the weighting "
     "function provided by the user, and n_j is a previously computed "
     "local density also provided by the user.");

  STRINGPC
      (symbol, m_symbolName,
       "Name of the symbol for the calculated local volume.");

  STRINGPC
      (densitySymbol, m_rhoSymbol,
       "Name of the local density to be used.");

  STRINGPC
      (weightingFunction, m_wfName,
       "Name of the weighting function to be used.");

  m_symbolName = "V";
  m_wfName = "default";
  m_rhoSymbol = "rho";

#ifdef _OPENMP
  m_particleCalculator = true;
#endif
}


void ValCalculatorVolume::setup()
{
  m_wf = M_SIMULATION->findWeightingFunction(m_wfName);

  m_cutoff = m_wf -> cutoff();

//   pair<bool, bool> persist(false, false);

  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("ValCalculatorVolume::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");

  if(m_allPairs)
  {

    size_t colour;
    // is the symbol already existing somewhere?
    for (colour = 0; colour < M_MANAGER->nColours(); ++colour)
    {
      if(Particle::s_tag_format[colour].attrExists(m_symbolName))
        throw gError("ValCalculatorVolume::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(colour) + "'. Second definition is not allowed with this Symbol calculator");
      if(!Particle::s_tag_format[colour].attrExists(m_rhoSymbol))
        throw gError("ValCalculatorVolume::setup", "Symbol " + m_rhoSymbol + " not found for species '" + M_MANAGER->species(colour) + "'.");
    }

        colour = 0;

//     if (!M_CONTROLLER->findIntegrator("IntegratorPosition", M_MANAGER->cp(0, 0)->firstSpecies()))
//       persist.first = true;

    // see CONVENTION5 for rule about persistency
    m_slots.first = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
    m_slots.second = m_slots.first;

    m_density_offset.first =
        Particle::s_tag_format[colour].offsetByName/*indexOf*/(m_rhoSymbol);
    m_density_offset.second =
        Particle::s_tag_format[colour].offsetByName/*indexOf*/(m_rhoSymbol);

    if(m_phaseUser == 0)
      Particle::registerCache_0
          (new ParticleCacheVolumeSelfContribution
          (colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
    if(m_phaseUser == 1)
      Particle::registerCache
          (new ParticleCacheVolumeSelfContribution
          (colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
    if(m_phaseUser == 2)
    {
      Particle::registerCache
          (new ParticleCacheVolumeSelfContribution
          (colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
      Particle::registerCache_0
          (new ParticleCacheVolumeSelfContribution
          (colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
    }

    FOR_EACH_COLOUR_PAIR
        (
        M_MANAGER,

//     bool newColour = false;

    m_density_offset.first =
        Particle::s_tag_format[cp->firstColour()].offsetByName/*indexOf*/(m_rhoSymbol);
    m_density_offset.second =
        Particle::s_tag_format[cp->secondColour()].offsetByName/*indexOf*/(m_rhoSymbol);

/*    persist.first = false;
    persist.second = false;*/


    // the following is old because replaced by CONVENTION5
#if 0
        // new rules
        // FIXME: not checked how meaningful this is for the case that both species don't have an IntegratorPosition
        // by the way, how meaningful is this case itself ?!?
    if (!M_CONTROLLER->findIntegrator("IntegratorPosition", cp->firstSpecies()) &&
         !M_CONTROLLER->findIntegrator("IntegratorPosition", cp->secondSpecies()))
    {
      persist.first = true;
      persist.second = true;
    }
#endif

    cp->setCutoff(m_cutoff);
    cp->setNeedPairs(true);

    if(colour < cp->firstColour() || colour < cp->secondColour())
    {
      ++colour;
      if(cp->firstColour() == colour)
      {
        // see CONVENTION5 for rule about persistency
        m_slots.first = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
        m_density_offset.first =
            Particle::s_tag_format[colour].offsetByName/*indexOf*/(m_rhoSymbol);

        if(m_phaseUser == 0)
          Particle::registerCache_0
              (new ParticleCacheVolumeSelfContribution(colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
        if(m_phaseUser == 1)
          Particle::registerCache
              (new ParticleCacheVolumeSelfContribution(colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
        if(m_phaseUser == 2)
        {
          Particle::registerCache_0
              (new ParticleCacheVolumeSelfContribution(colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
          Particle::registerCache
              (new ParticleCacheVolumeSelfContribution(colour, m_slots.first, m_symbolName, m_density_offset.first, m_wf));
        }
      }
      if(cp->secondColour() == colour)
      {
        // see CONVENTION5 for rule about persistency
        m_slots.second = Particle::s_tag_format[colour].addAttribute(m_symbolName, m_datatype, /*persist.second*/false, m_symbolName).offset;
        m_density_offset.second =
            Particle::s_tag_format[colour].offsetByName/*indexOf*/(m_rhoSymbol);
        if(m_phaseUser == 0)
          Particle::registerCache_0
              (new ParticleCacheVolumeSelfContribution
              (colour, m_slots.second, m_symbolName, m_density_offset.second, m_wf));
        if(m_phaseUser == 1)
          Particle::registerCache
              (new ParticleCacheVolumeSelfContribution
              (colour, m_slots.second, m_symbolName, m_density_offset.second, m_wf));
        if(m_phaseUser == 2)
        {
          Particle::registerCache
              (new ParticleCacheVolumeSelfContribution
              (colour, m_slots.second, m_symbolName, m_density_offset.second, m_wf));
          Particle::registerCache_0
              (new ParticleCacheVolumeSelfContribution
              (colour, m_slots.second, m_symbolName, m_density_offset.second, m_wf));
        }
      }
    }     // end: newColour = true;
    else
    {
      // here we have to make sure at least that the offsets are correct
      m_slots.first = Particle::s_tag_format[cp->firstColour()].offsetByName(m_symbolName);
      m_density_offset.first =
          Particle::s_tag_format[cp->firstColour()].offsetByName(m_rhoSymbol);
      m_slots.second = Particle::s_tag_format[cp->secondColour()].offsetByName(m_symbolName);
      m_density_offset.second =
          Particle::s_tag_format[cp->secondColour()].offsetByName(m_rhoSymbol);
    }



    vector<ColourPair*>::iterator cpTester = __cp;
        // is it the last calculator to be created?
    if(++cpTester == __end)
    {
      if(m_phaseUser == 0)
        cp->registerCalc_0(this);
      else if(m_phaseUser == 0)
        cp->registerCalc(this);
      else // so it is 2
      {
        ValCalculator* vc = /*copyMySelf()*/new ValCalculatorVolume(*this);
        assert(((ValCalculatorVolume*) (vc))->m_symbolName == m_symbolName);
        assert(vc->stage() == m_stage);
        assert(((ValCalculatorVolume*) (vc))->m_wf == m_wf);
        assert(((ValCalculatorVolume*) (vc))->m_wfName == m_wfName);
        assert(((ValCalculatorVolume*) (vc))->m_rhoSymbol == m_rhoSymbol);
        assert(((ValCalculatorVolume*) (vc))->m_slots.first == m_slots.first);
        assert(((ValCalculatorVolume*) (vc))->m_slots.second == m_slots.second);
        assert(((ValCalculatorVolume*) (vc))->m_parent == m_parent);
        assert(((ValCalculatorVolume*) (vc))->m_species.first == m_species.first);
        assert(((ValCalculatorVolume*) (vc))->m_species.second == m_species.second);

        cp->registerCalc(vc);
        cp->registerCalc(this);

      }
    }
        // No? Then make a copy
    else
    {
      ValCalculator* vc = /*copyMySelf()*/new ValCalculatorVolume(*this);
      assert(((ValCalculatorVolume*) (vc))->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((ValCalculatorVolume*) (vc))->m_wf == m_wf);
      assert(((ValCalculatorVolume*) (vc))->m_wfName == m_wfName);
      assert(((ValCalculatorVolume*) (vc))->m_rhoSymbol == m_rhoSymbol);
      assert(((ValCalculatorVolume*) (vc))->m_slots.first == m_slots.first);
      assert(((ValCalculatorVolume*) (vc))->m_slots.second == m_slots.second);
      assert(((ValCalculatorVolume*) (vc))->m_parent == m_parent);
      assert(((ValCalculatorVolume*) (vc))->m_species.first == m_species.first);
      assert(((ValCalculatorVolume*) (vc))->m_species.second == m_species.second);

      if(m_phaseUser == 0)
        cp->registerCalc_0(vc);
      else if(m_phaseUser == 0)
        cp->registerCalc(vc);
      else // so it is 2
      {
        cp->registerCalc(vc);

        vc = /*copyMySelf()*/new ValCalculatorVolume(*this);
        assert(((ValCalculatorVolume*) (vc))->m_symbolName == m_symbolName);
        assert(vc->stage() == m_stage);
        assert(((ValCalculatorVolume*) (vc))->m_wf == m_wf);
        assert(((ValCalculatorVolume*) (vc))->m_wfName == m_wfName);
        assert(((ValCalculatorVolume*) (vc))->m_rhoSymbol == m_rhoSymbol);
        assert(((ValCalculatorVolume*) (vc))->m_slots.first == m_slots.first);
        assert(((ValCalculatorVolume*) (vc))->m_slots.second == m_slots.second);
        assert(((ValCalculatorVolume*) (vc))->m_parent == m_parent);
        assert(((ValCalculatorVolume*) (vc))->m_species.first == m_species.first);
        assert(((ValCalculatorVolume*) (vc))->m_species.second == m_species.second);

        cp->registerCalc_0(vc);
      }
    }
        );
  }
  else /*m_allPairs == false*/
  {
    if(m_species.first == "undefined")
      throw gError("ValCalculatorVolume::setup", "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");
    if(m_species.second == "undefined")
      throw gError("ValCalculatorVolume::setup", "Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled.");

    ColourPair* cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

    cp->setCutoff(m_cutoff);
    cp->setNeedPairs(true);

    if(Particle::s_tag_format[cp->firstColour()].attrExists(m_symbolName))
      throw gError("ValCalculatorVolume::setup", "Symbol " + m_symbolName + " already existing for species '" + M_MANAGER->species(cp->firstColour()) + "'. Second definition is not allowed for this Calculator.");
    if(Particle::s_tag_format[cp->secondColour()].attrExists(m_symbolName))
      throw gError("ValCalculatorVolume::setup", "Symbol " + m_symbolName + " already existing for species '" + M_MANAGER->species(cp->secondColour()) + "'. Second definition is not allowed for this Calculator.");
    if(!Particle::s_tag_format[cp->firstColour()].attrExists(m_rhoSymbol))
      throw gError("ValCalculatorVolume::setup", "Symbol " + m_rhoSymbol + " not found for species '" + M_MANAGER->species(cp->firstColour()) + "'.");
    if(!Particle::s_tag_format[cp->secondColour()].attrExists(m_rhoSymbol))
      throw gError("ValCalculatorVolume::setup", "Symbol " + m_rhoSymbol + " not found for species '" + M_MANAGER->species(cp->secondColour()) + "'.");

    // the following is old because of CONVENTION5
#if 0
    // new rules
    // FIXME: not checked how meaningful this is for the case that both species don't have an IntegratorPosition
    // by the way, how meaningful is this case itself ?!?
    if (!M_CONTROLLER->findIntegrator("IntegratorPosition", cp->firstSpecies()) &&
         !M_CONTROLLER->findIntegrator("IntegratorPosition", cp->secondSpecies()))
    {
      persist.first = true;
      persist.second = true;
    }
#endif
    // see CONVENTION5 for rule about persistency
    m_slots.first = Particle::s_tag_format[cp->firstColour()].addAttribute(m_symbolName, m_datatype, /*persist.first*/false, m_symbolName).offset;
    m_density_offset.first =
        Particle::s_tag_format[cp->firstColour()].offsetByName/*indexOf*/(m_rhoSymbol);

    if(m_phaseUser == 0)
      Particle::registerCache_0
          (new ParticleCacheVolumeSelfContribution
          (cp->firstColour(), m_slots.first, m_symbolName, m_density_offset.first, m_wf));
    if(m_phaseUser == 1)
      Particle::registerCache
          (new ParticleCacheVolumeSelfContribution
          (cp->firstColour(), m_slots.first, m_symbolName, m_density_offset.first, m_wf));
    if(m_phaseUser == 2)
    {
      Particle::registerCache_0
          (new ParticleCacheVolumeSelfContribution
          (cp->firstColour(), m_slots.first, m_symbolName, m_density_offset.first, m_wf));
      Particle::registerCache
          (new ParticleCacheVolumeSelfContribution
          (cp->firstColour(), m_slots.first, m_symbolName, m_density_offset.first, m_wf));
    }

    if(cp->firstColour() != cp->secondColour())
    {
      // see CONVENTION5 for rule about persistency
      m_slots.second = Particle::s_tag_format[cp->secondColour()].addAttribute(m_symbolName, m_datatype, /*persist.second*/false, m_symbolName).offset;
      m_density_offset.second =
          Particle::s_tag_format[cp->secondColour()].offsetByName/*indexOf*/(m_rhoSymbol);

      if(m_phaseUser == 0)
        Particle::registerCache_0
            (new ParticleCacheVolumeSelfContribution
            (cp->secondColour(), m_slots.second, m_symbolName, m_density_offset.second, m_wf));
      if(m_phaseUser == 1)
        Particle::registerCache
            (new ParticleCacheVolumeSelfContribution
            (cp->secondColour(), m_slots.second, m_symbolName, m_density_offset.second, m_wf));
      if(m_phaseUser == 2)
      {
        Particle::registerCache_0
            (new ParticleCacheVolumeSelfContribution
            (cp->secondColour(), m_slots.second, m_symbolName, m_density_offset.second, m_wf));
        Particle::registerCache
            (new ParticleCacheVolumeSelfContribution
            (cp->secondColour(), m_slots.second, m_symbolName, m_density_offset.second, m_wf));
      }
    }
    else m_slots.second = m_slots.first;

    if(m_phaseUser == 0)
      cp->registerCalc_0(this);
    else if(m_phaseUser == 1)
      cp->registerCalc(this);
    else // so it is 2
    {
      ValCalculator* vc = /*copyMySelf()*/new ValCalculatorVolume(*this);
      assert(((ValCalculatorVolume*) (vc))->m_symbolName == m_symbolName);
      assert(vc->stage() == m_stage);
      assert(((ValCalculatorVolume*) (vc))->m_wf == m_wf);
      assert(((ValCalculatorVolume*) (vc))->m_wfName == m_wfName);
      assert(((ValCalculatorVolume*) (vc))->m_rhoSymbol == m_rhoSymbol);
      assert(((ValCalculatorVolume*) (vc))->m_slots.first == m_slots.first);
      assert(((ValCalculatorVolume*) (vc))->m_slots.second == m_slots.second);
      assert(((ValCalculatorVolume*) (vc))->m_parent == m_parent);
      assert(((ValCalculatorVolume*) (vc))->m_species.first == m_species.first);
      assert(((ValCalculatorVolume*) (vc))->m_species.second == m_species.second);

      cp->registerCalc(this);
      cp->registerCalc_0(vc);
    }

  }
}


#ifdef _OPENMP
void ValCalculatorVolume::mergeCopies(ColourPair* cp, int thread_no) {
  size_t slot1 = m_slots.first;
  size_t slot2 = m_slots.second;

  size_t copySlot1 = m_copy_slots[thread_no].first;
  size_t copySlot2 = m_copy_slots[thread_no].second;
  size_t vecSlot1 = m_vector_slots.first;
  size_t vecSlot2 = m_vector_slots.second;

  FOR_EACH_PARTICLE_C
   (M_PHASE, cp->firstColour(),
      __iSLFE->tag.doubleByOffset(slot1) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1];
        (*__iSLFE->tag.vectorDoubleByOffset(copySlot1))[vecSlot1] = 0;
   );
  FOR_EACH_PARTICLE_C
   (M_PHASE, cp->secondColour(),
      __iSLFE->tag.doubleByOffset(slot2) += (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2];
        (*__iSLFE->tag.vectorDoubleByOffset(copySlot2))[vecSlot2] = 0;
   );
}
#endif




