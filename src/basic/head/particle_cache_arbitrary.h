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



#ifndef __PARTICLE_CACHE_ARBITRARY_H
#define __PARTICLE_CACHE_ARBITRARY_H

#include "particle_cache.h"

#include "function_particle.h"
#include "symbol.h"

/*!
 * Functions to cache completely user-defined properties for the 
 * particles, which depend only on other particle properties.
 */
class ParticleCacheArbitrary : public ParticleCache
{
  protected:
  
    /*!
    * The mathematical expression to be computed
    */
    string m_expression;
    
    /*!
    * The \a FunctionParticle computing \a m_expression
    */
    FunctionParticle* m_function;
    
    /*!
     * Initialise the PropertyList.
     */
    void init();
  
    /*!
     * Helper for setting m_offset
     */
    virtual void setupOffset();
    
    /*!
     * Return a copy of the current object
     */
    virtual ParticleCache* copyMySelf() = 0;
     
    /*!
     * Adds the expressions used by this \a Symbol to the given list. 
     * @param usedSymbols List to be filled with own instances of \a TypedValue
     */
    virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols);

    /*!
     * Returns the strings of those \a Symbols that the given class depends on
     * due to hard-coded reasons (not due to runtime compiled expressions).
     * @param usedSymbols List to add the strings to.
     */
    virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
    {

    }

    /*!
     * Helper function for setting the return type of \a m_function
     */
    virtual void setFunctionReturnType() = 0;
    
  public:
    /*!
     * Constructor
     */
    ParticleCacheArbitrary(/*Node*/Simulation* parent);

    /*!
     * Destructor
     */
    virtual ~ParticleCacheArbitrary();

    /*!
     * Setup this ParticleCache
     */
    virtual void setup();

    /*!
     * Delete old function and set \a m_function to the new one given
     * @param[in] Pointer to new function
     */
    virtual void setFunctionParticle(FunctionParticle* fp) {

      delete m_function;

      m_function = fp;

      setFunctionReturnType();
    }
};

#endif
