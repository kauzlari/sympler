/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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



#ifndef __PARTICLE_CACHE_H
#define __PARTICLE_CACHE_H

#include "particle.h"

#include "symbol.h"


/*!
 * Classes to cache properties for the particles. Similar to the ValCalculators for
 * the pairs.
 */
class ParticleCache : public Symbol
{

 protected:
  
  /*!
   * memory offset where the computed attribute is located
   */
  size_t m_offset;
  
  /*!
   * Color of the particles this cache is used for.
   */
  string m_species;
  
  /*!
   * Color of the particles this cache is used for.
   */
  size_t m_colour;
  
  /*!
   * Initialise the PropertyList.
   */
  void init();

  /*!
   * If it belongs to a Node structure, setup this instance of
   * \a ParticleCache
   * FIXME: This function was added pretty recently (2018-02-19), and 
   * most of the children do not yet call it in their own setup()
   */
  virtual void setup();

  /*!
   * Checks existence of output symbol given by \a m_symbolName and 
   * reacts according to values for \a m_overwrite, etc.
   * @param colour Particle colour to be checked
   */
  virtual void checkOutputSymbolExistence(size_t colour);

  /*!
   * Checks existence of input symbols required by this \a ParticleCache
   * in a hard-coded fashion (i.e., not through runtime compiled 
   * expressions).
   * Default behaviour defined here: no input symbols
   * @param colour Particle colour to be checked
   */
  virtual void checkInputSymbolExistences(size_t colour) {
  }

  /*!
   * Helper function which removes brackets from single terms in \a Function. 
   * expressions. E.g.: "[r]" becomes "r"
   * Differently to class \a ValCalculator, there shouldn't be any indices to 
   * be removed
   * @param name Single term from a \a Function expression
   */
  virtual void cleanSymbol(string& name) const;

  /*!
   * Return a copy of the current object
   * Should be pure virtual if class is abstract anyways.
   */
  virtual ParticleCache* copyMySelf() = 0;

  
 public:
  /*!
   * Constructor for Node hierarchy
   */
  ParticleCache(/*Node*/Simulation* parent);
  
  /*!
   * Constructor
   */
  ParticleCache(size_t colour, size_t offset, string symbolName);
  
  /*!
   * Destructor
   */
  virtual ~ParticleCache();

  /*!
   * Compute the cache for particle \a p
   * @param p The particle to compute values for
   */
  virtual void computeCacheFor(Particle* p) = 0;

  /*!
   * Register the additional degrees of freedom with the \a Particle s \a DataFormat
   */
  virtual void registerWithParticle() = 0;

  /*!
   * Are those two caches identical?
   * @param c Is this the same cache?
   */
  virtual bool operator==(const ParticleCache &c) const = 0;

  /*!
   * Return the color this cache is for
   */
  size_t colour() const {
    return m_colour;
  }

  /*!
   * Return the color this cache writes into. HEre, the default behaviour is implemented
   */
  virtual size_t writeColour() const {
    return m_colour;
  }

  /*!
   * returns the offset
   */
  size_t offset() const {
    return m_offset; 
  }
  
};

#endif
