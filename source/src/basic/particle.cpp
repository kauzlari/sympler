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



#include "particle.h"

#include "misc.h"
#include "particle_cache.h"

using namespace std;

/*--- Particle ---*/

vector<Data/*Format*/> Particle::s_tag_format;
vector<vector<vector<ParticleCache*> > > Particle::s_cached_properties;
vector<ParticleCache*> Particle::s_cached_flat_properties;
size_t Particle::s_maxStage = 0;
vector<vector<vector<ParticleCache*> > > Particle::s_cached_properties_0;
vector<ParticleCache*> Particle::s_cached_flat_properties_0;
size_t Particle::s_maxStage_0 = 0;

ParticleCache *Particle::registerCache(ParticleCache *c)
{
  ParticleCache *cc = NULL;
  
//   size_t stage = c->stage();

  FOR_EACH
    (vector<ParticleCache*>,
     s_cached_flat_properties,
     if ((*(*__iFE)) == (*c))
       cc = *__iFE;
     );

  // FIXME: Next could cause trouble at some point
  if (cc)
    throw gError("Particle::registerCache", "Cache " + c->name() + " already registered. Contact the programmer.");
    //     return cc;

  c->registerWithParticle();
  
//   s_cached_properties[c->colour()][stage].push_back(c);
  s_cached_flat_properties.push_back(c);

  MSG_DEBUG("Particle::registerCache", "registering cache: " << c->className() << ", symbol=" << c->mySymbolName());

  return c;
}

ParticleCache *Particle::registerCache_0(ParticleCache *c)
{
  ParticleCache *cc = NULL;
  
  FOR_EACH
    (vector<ParticleCache*>,
     s_cached_flat_properties_0,
     if ((*(*__iFE)) == (*c)) cc = *__iFE;
     );

  // FIXME: Next could cause trouble at some point
  if (cc)
    throw gError("Particle::registerCache", "Cache " + c->name() + " already registered. Contact the programmer.");

  c->registerWithParticle();
  
  s_cached_flat_properties_0.push_back(c);

  return c;
}

void Particle::sortStages()
{
  FOR_EACH
    (vector<ParticleCache*>,
     s_cached_flat_properties,
     int stage = (*__iFE)->stage();
     assert(stage > -1);
     if((size_t) stage > s_maxStage) s_maxStage = stage;
     if(s_cached_properties[(*__iFE)->colour()].size() <= (size_t) stage)
       s_cached_properties[(*__iFE)->colour()].resize(stage+1);
     s_cached_properties[(*__iFE)->colour()][stage].push_back(*__iFE);
     );
  MSG_DEBUG("Particle::sortStages", "highest stage: " << s_maxStage);  
}

void Particle::sortStages_0()
{
  FOR_EACH
    (vector<ParticleCache*>,
     s_cached_flat_properties_0,
     int stage = (*__iFE)->stage();
     assert(stage > -1);
     if((size_t) stage > s_maxStage_0) s_maxStage_0 = stage;
     if(s_cached_properties_0[(*__iFE)->colour()].size() <= (size_t) stage)
       s_cached_properties_0[(*__iFE)->colour()].resize(stage+1);
     s_cached_properties_0[(*__iFE)->colour()][stage].push_back(*__iFE);
     );
  MSG_DEBUG("Particle::sortStages_0", "highest stage: " << s_maxStage_0);  
  
}
