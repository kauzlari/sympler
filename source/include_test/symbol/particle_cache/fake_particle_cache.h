/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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



#ifndef __FAKE_PARTICLE_CACHE_H
#define __FAKE_PARTICLE_CACHE_H

#include "particle_cache.h"

/*!
 * Fake class derived from class \a FakeParticleCache for testing purposes. this
 * class mainly implements default behaviour for the purely virtual functions
 */
class FakeParticleCache : public ParticleCache
{

 protected:
  
  /*!
   * Return a copy of the current object.
   * Here, default fake behaviour is implemented which will NOT be tested.
   * Function was pure virtual in parent class \a ParticleCache .
   */
  virtual FakeParticleCache* copyMySelf() {
  	return this;
  }

  
 public:
  /*!
   * Constructor for Node hierarchy
   */
  FakeParticleCache(/*Node*/Simulation* parent);
  
  /*!
   * Constructor
   */
  FakeParticleCache(size_t colour, size_t offset, string symbolName);
  
  /*!
   * Destructor
   */
  virtual ~FakeParticleCache();

  /*!
   * Calls the protected \a ParticleCache::cleanSymbol() for test purposes.
   * @param[in] name Symbol name to be cleaned.
   */
  virtual void cleanSymbolPublic(string& name) const {

  	cleanSymbol(name);
  }

  /*!
   * Sets protected member \a ParticleCache::m_datatype for testing purposes
   */
  virtual void setDatatype(DataFormat::datatype_t datatype) {

  	m_datatype = datatype;
  }

  /*!
   * Sets protected member \a ParticleCache::m_overwrite for testing purposes
   */
  virtual void setOverwrite(bool overwrite) {

  	m_overwrite = overwrite;
  }

  /*!
   * Compute the cache for particle \a p
   * Here, default fake behaviour is implemented which will NOT be tested.
   * Function was pure virtual in parent class \a ParticleCache .
   * @param p The particle to compute values for
   */
  virtual void computeCacheFor(Particle* p) {}

  /*!
   * Register the additional degrees of freedom with the \a Particle 's
   * \a DataFormat
   * Here, default fake behaviour is implemented which will NOT be tested.
   * Function was pure virtual in parent class \a ParticleCache .
   */
  virtual void registerWithParticle() {}

  /*!
   * Are those two caches identical?
   * Here, default fake behaviour is implemented which will NOT be tested.
   * Function was pure virtual in parent class \a ParticleCache .
   * @param c Is this the same cache?
   */
  virtual bool operator==(const ParticleCache &c) const {

  	return false;

  }

};

#endif
