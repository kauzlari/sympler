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


#include "transfer_particle_vector.h"

#include "simulation.h"
#include "manager_cell.h"

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()

const SymbolRegister<TransferParticleVector> transfer_particle_vector("TransferParticleVector");


TransferParticleVector::TransferParticleVector(Simulation* parent)
  : ParticleCacheArbitrary(parent), m_pointerToParticleList(NULL), m_targetColour(100000)
{
  m_function.setReturnType(Variant::VECTOR);
  m_datatype = DataFormat::POINT;
  init();
}

TransferParticleVector::~TransferParticleVector()
{
}

void TransferParticleVector::init()
{
  m_properties.setClassName("TransferParticleVector");
  m_properties.setDescription("User defined vector-type particle property of a 'targetSpecies', which depends exclusively on properties of particles of a \"source species\" given by the attribute 'species'. It may be defined by the attribute 'expression'.");

  STRINGPC
      (targetSpecies, m_targetSpecies,
       "Name of the target species.");

  m_targetSpecies = "undefined";

  STRINGPC 
      (targetSlot, m_targetSlotName,
       "Name for the symbol of the \"source particle\" storing the slot of the \"target particle\" it delivers data to.");

  m_targetSlotName = "undefined";

}

void TransferParticleVector::setup()
{
   MSG_DEBUG("TransferParticleVector::setup", "START");


  M_SIMULATION->controller()->registerForSetupAfterParticleCreation(this);

  if(m_targetSpecies == "undefined")
    throw gError("TransferParticleVector::setup", "Attribute 'targetSpecies' has value \"undefined\""); 
  if(m_targetSlotName == "undefined")
    throw gError("TransferParticleVector::setup", "'targetSlot' has value \"undefined\"");

  m_targetColour = M_PHASE->manager()->getColour(m_targetSpecies);

  ParticleCacheArbitrary::setup();
   MSG_DEBUG("TransferParticleVector::setup", "END");
}

// called from ParticleCacheArbitrary::setup()
void TransferParticleVector::setupOffset()
{
  MSG_DEBUG("TransferParticleVector::setupOffset", "START");
  if(m_overwrite)
  { 
    // so the attribute should already exist
    try
    {
      m_offset = Particle::s_tag_format[m_targetColour].indexOf(m_symbolName, m_datatype);
      m_offset = Particle::s_tag_format[m_targetColour].offsetByIndex(m_offset);
    }
    catch(gError& err)
    {
      throw gError("ParticleCache::setup", "search for symbol failed. The message was " + err.message()); 
    }
  }
  else
  {
    // so the attribute shouldn't yet exist
    if(Particle::s_tag_format[m_targetColour].attrExists(m_symbolName))
      throw gError("ParticleCache::setup", "Symbol " + m_symbolName + " already existing. Second definition is not allowed for 'overwrite = \"no\"'.");
    else 
      // FIXME: let's try if it works for all cases that persistency = false
      m_offset = Particle::s_tag_format[m_targetColour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
    MSG_DEBUG("ParticleCacheArbitrary::setupOffset", "Offset for " << m_symbolName << " = " << m_offset);
  } 

  // if the transfer-slot symbol doesn't yet exist, we add it
  if(Particle::s_tag_format[m_colour].attrExists(m_targetSlotName)) {
    m_targetSlotOffset = Particle::s_tag_format[m_colour].indexOf(m_targetSlotName, DataFormat::INT);
    m_targetSlotOffset = Particle::s_tag_format[m_colour].offsetByIndex(m_targetSlotOffset);
  }
  else {
    // This MUST be persistent !!!
    // Currently (2009-09-06) we create the symbol here, and then, no module writes it except ParticleCreatorFile once at the beginning
    m_targetSlotOffset = Particle::s_tag_format[m_colour].addAttribute(m_targetSlotName, DataFormat::INT, true, m_targetSlotName).offset;
  MSG_DEBUG("TransferParticleVector::setupOffset", "Offset for " << m_targetSlotName << " = " << m_targetSlotOffset);
  }

  MSG_DEBUG("TransferParticleVector::setupOffset", "END");

}

void TransferParticleVector::setupAfterParticleCreation()
{
  // set the helper
  m_pointerToParticleList = &(M_PHASE->particles(m_targetColour));

  // make a consistency check, whether all required target particles really exist
  FOR_EACH_FREE_PARTICLE_C
    (
     M_PHASE,
     m_colour,
     size_t targetSlot = __iSLFE->tag.intByOffset(m_targetSlotOffset);
     if(targetSlot >= m_pointerToParticleList->size())
       throw gError("TransferParticleVector::setupAfterParticleCreation", "ERROR: Source particle at position (" + ObjToString(__iSLFE->r.x) + ", " + ObjToString(__iSLFE->r.y) + ", " + ObjToString(__iSLFE->r.z) + ") requires non-existent target particle with slot " + ObjToString(targetSlot)+ ". There are just " + ObjToString(m_pointerToParticleList->size()) + " target particles which have been created. Check this!" );
     );
}
