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


#ifndef __PARTICLE_RAND_NORM_VECTOR_H
#define __PARTICLE_RAND_NORM_VECTOR_H

#include "particle_cache_arb_rng.h"
#include "random.h"

/*!
 * User-defined symbol for a particle computing a vector of 3 statistically independent random normaly distributed scalars multiplied by 
 * previously computed particle properties (symbols).
 */
class ParticleRandNormVector : public ParticleCacheArbRNG
{
  protected:
     
    /*!
     * Initialise the PropertyList.
     */
    void init();
  
    
    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new ParticleRandNormVector(*this);
    }
    
  public:
  /*!
   * Constructor
   */
    ParticleRandNormVector(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ParticleRandNormVector();

    /*!
     * Precompute the random numbers for each particle with zero mean and unit variance
     */
    virtual void precompute()
    {
      point_t average = {{{0,0,0}}};
      point_t sumOfSquares = {{{0,0,0}}};
      point_t sumOfSquaresCentred = {{{0,0,0}}};
      point_t* thisPoint;

      // draw random numbers
      FOR_EACH_FREE_PARTICLE_C
	(m_phasePointer, m_colour,
	 thisPoint = &(__iSLFE->tag.pointByOffset(m_offset));
	 for(size_t i = 0; i < SPACE_DIMS; ++i) {
	   (*thisPoint)[i] = m_rng.normal(1);
	   average[i] += (*thisPoint)[i];
	   sumOfSquares[i] += ((*thisPoint)[i])*((*thisPoint)[i]);
	 }
	 );

      // shifting the average to zero and rescaling the variance 
      // is only done if the number of particles is at least 
      // the user-specified limit (default=2)
      size_t nOfP = m_phasePointer->returnNofPartC(m_colour);

      // size_t cast is safe because m_plimit was checked before
      if(nOfP >= size_t(m_plimit)) {
	// compute average
	for(size_t i = 0; i < SPACE_DIMS; ++i) {
	  average[i] /= nOfP;
	}

	// subtract average
	FOR_EACH_FREE_PARTICLE_C
	  (m_phasePointer, m_colour,
	   thisPoint = &(__iSLFE->tag.pointByOffset(m_offset));
	   for(size_t i = 0; i < SPACE_DIMS; ++i) {
	     (*thisPoint)[i] -= average[i];
	     sumOfSquaresCentred[i] += ((*thisPoint)[i])*((*thisPoint)[i]);
	   }
	   );
	
	point_t scaling;
	for(size_t i = 0; i < SPACE_DIMS; ++i) {
	  scaling[i] = sqrt(sumOfSquares[i]/sumOfSquaresCentred[i]);
	}

	// rescale to original variance
	FOR_EACH_FREE_PARTICLE_C
	  (m_phasePointer, m_colour,
	   thisPoint = &(__iSLFE->tag.pointByOffset(m_offset));
	   for(size_t i = 0; i < SPACE_DIMS; ++i) {
	     (*thisPoint)[i] *= scaling[i];
	   }
	   );	
      }

    }

      /*!
     * Compute the cache for particle \a p
     * @param p The particle to compute values for
       */
    virtual void computeCacheFor(Particle* p)
    {
      
      // new style: random number (zero mean, unit variance) was already precomputed in the corresponding tag-slot
      point_t* vec = &(p->tag.pointByOffset(m_offset));
/*       m_function(vec, p); */
      point_t temp;
      m_function(&temp, p);
      for(size_t i = 0; i < SPACE_DIMS; ++i)
	(*vec)[i] *= temp[i];
/* 	(*vec)[i] *= m_rng.normal(1); */

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
        ParticleRandNormVector *cc = (ParticleRandNormVector*) &c;

        return
          (
          m_expression == cc->m_expression &&
          m_overwrite == cc->m_overwrite &&
          m_stage == cc->m_stage &&
          m_colour == cc->m_colour &&
          m_symbolName == cc->m_symbolName &&
          m_datatype == cc->m_datatype
          );
      }
      else return false;
    }
    
  /*!
     * Setup this ParticleCache
   */
    virtual void setup();
};

#endif
