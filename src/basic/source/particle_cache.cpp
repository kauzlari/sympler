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



#include "particle_cache.h"

ParticleCache::ParticleCache(/*Node*/Simulation* parent/*size_t colour*/)
/*: m_colour(colour), m_stage(0)*/
  : Symbol(parent)
{
  init();
}

ParticleCache::ParticleCache(size_t colour, size_t offset, string symbolName)
  : Symbol(symbolName)
{
  m_colour = colour;
  m_offset = offset;
//   m_symbolName = symbolName;
}

ParticleCache::~ParticleCache()
{
}

void ParticleCache::init()
{
  m_properties.setClassName("ParticleCache");

  STRINGPC
      (species, m_species,
       "Name for the species of the particles, this Symbol is used for. If set to \"ALL\", the Symbol will be used for all registered species.");
  
  m_species = "undefined";
   
}
