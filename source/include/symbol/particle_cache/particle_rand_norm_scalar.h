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


#ifndef __PARTICLE_RAND_NORM_SCALAR_H
#define __PARTICLE_RAND_NORM_SCALAR_H

#include "particle_cache_arb_rng.h"
#include "random.h"

/*!
 * User-defined random normaly distributed scalar symbol for a particle multiplied by 
 * previously computed particle properties (symbols).
 */
class ParticleRandNormScalar : public ParticleCacheArbRNG
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
      return new ParticleRandNormScalar(*this);
    }

    /*!
     * Helper function for setting the return type of \a m_function
     */
    virtual void setFunctionReturnType(){
      m_function->setReturnType(Variant::SCALAR);
    }
    
  public:
  /*!
   * Constructor
   */
    ParticleRandNormScalar(/*Node*/Simulation* parent);

  /*!
     * Destructor
   */
    virtual ~ParticleRandNormScalar();

    /*!
     * Precompute the random numbers for each particle with zero mean and unit variance
     */
    virtual void precompute()
    {
      double average = 0;
      double sumOfSquares = 0;
      double sumOfSquaresCentred = 0;
      double* thisNum;

      // draw random numbers
      FOR_EACH_FREE_PARTICLE_C
	(m_phasePointer, m_colour,
	 thisNum = &(__iSLFE->tag.doubleByOffset(m_offset));
	 *thisNum = m_rng.normal(1);
	 average += *thisNum;
	 sumOfSquares += (*thisNum)*(*thisNum);
	 );

      // shifting the average to zero and rescaling the variance 
      // is only done if the number of particles is at least 
      // the user-specified limit (default=2)
      size_t nOfP = m_phasePointer->returnNofPartC(m_colour);

      // size_t cast is safe because m_plimit was checked before
      if(nOfP >= size_t(m_plimit)) {
	// compute average
	average /= nOfP;
	
	// subtract average
	FOR_EACH_FREE_PARTICLE_C
	  (m_phasePointer, m_colour,
	   thisNum = &(__iSLFE->tag.doubleByOffset(m_offset));
	   *thisNum -= average;
	   sumOfSquaresCentred += (*thisNum)*(*thisNum);
	   );
	
	double scaling = sqrt(sumOfSquares/sumOfSquaresCentred);
	
	// rescale to original variance
	FOR_EACH_FREE_PARTICLE_C
	  (m_phasePointer, m_colour,
	   __iSLFE->tag.doubleByOffset(m_offset) *= scaling;
	   );
	
      }
    }


      /*!
     * Compute the cache for particle \a p
     * @param p The particle to compute values for
       */
    virtual void computeCacheFor(Particle* p)
    {
      
      double temp;
      (*m_function)(&temp, p);

      // new style: random number (zero mean, unit variance) was already precomputed in the corresponding tag-slot
      p->tag.doubleByOffset(m_offset) *= temp;
/*       p->tag.doubleByOffset(m_offset) = m_rng.normal(1)*temp; */
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
        ParticleRandNormScalar *cc = (ParticleRandNormScalar*) &c;

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
