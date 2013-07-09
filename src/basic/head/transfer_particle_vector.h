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


#ifndef __TRANSFER_PARTICLE_VECTOR_H
#define __TRANSFER_PARTICLE_VECTOR_H

#include "particle_cache_arbitrary.h"

#include "particle_list.h"

/*!
 * User-defined vector symbol for a particle of a "target-species" 
 * that depends only on previously computed particle properties
 * of a "source-species".
 */
class TransferParticleVector : public ParticleCacheArbitrary
{
  protected:

  /*!
   * offset to the "source-particle's" symbol representing the 
   * slot of the corresponding "target-particle" 
   */
  size_t m_targetSlotOffset;

  /*!
   * name of the "target" species
   * slot of the corresponding "target-particle" 
   */
  string m_targetSpecies;

  /*!
   * name of the symbol representing the
   */
  string m_targetSlotName;
  
  /*!
   * helper pointing to the right "target"-particle list
   */
  ParticleList* m_pointerToParticleList;

  /*!
   * helper for storing the "target"-colour corresponding to \a m_targetSpecies
   */
  size_t m_targetColour;

    /*!
   * Initialise the PropertyList.
     */
    void init();
      
    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new TransferParticleVector(*this);
    }

    /*!
     * Setup the offsets calling ParticleCacheArbitrary::setupOffset() and the additional ones of this class
     */
    virtual void setupOffset();
    
  public:
  /*!
   * Constructor
   */
    TransferParticleVector(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~TransferParticleVector();

      /*!
     * Compute the cache for particle \a p
     * @param p The particle to compute values for
       */
    virtual void computeCacheFor(Particle* sourceP)
    { 
      size_t targetSlot = sourceP->tag.intByOffset(m_targetSlotOffset);

      Particle* targetP = &((*m_pointerToParticleList)[targetSlot]);

      point_t temp;

      m_function(&temp, sourceP);

/*       if(targetP->mySlot == 100) */
/* 	MSG_DEBUG("TransferParticleVector::computeCacheFor", "p100:BEFORE " << m_expression << ", " + m_symbolName << " = " << targetP->tag.pointByOffset(m_offset) << "stage=" << m_stage); */
      
      targetP->tag.pointByOffset(m_offset) += temp;

/*       if(targetP->mySlot == 100) */
/* 	MSG_DEBUG("TransferParticleVector::computeCacheFor", "p100:AFTER " << m_expression << ", " + m_symbolName << " = " << targetP->tag.pointByOffset(m_offset) << "stage=" << m_stage); */
      
    }

  /*!
     * Register the additional degrees of freedom with the \a Particle s \a DataFormat
   */
    virtual void registerWithParticle()
    {
    }
    
  /*!
     * Are those two caches identical?
     * @param c Is this the same cache?
   */
    virtual bool operator==(const ParticleCache &c) const
    {
      if (typeid(c) == typeid(*this)) 
      {
        TransferParticleVector *cc = (TransferParticleVector*) &c;

        return
            (
            m_expression == cc->m_expression &&
            m_overwrite == cc->m_overwrite &&
            m_stage == cc->m_stage &&
            m_colour == cc->m_colour &&
            m_symbolName == cc->m_symbolName &&
            m_datatype == cc->m_datatype &&
            m_offset == cc->m_offset
            );
      }
      else return false;
    }
    
  /*!
     * Setup this ParticleCache
   */
    virtual void setup();

    /*!
     * Additional setup
     */
    virtual void setupAfterParticleCreation();

    /*!
     * Return the colour this cache writes into. Here, the default behaviour is overwritten
     */
    virtual size_t writeColour() const {
      return m_targetColour;
    }


};

#endif
