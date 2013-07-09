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



#include "particle_list.h"



/*--- ParticleList ---*/

//ParticleList::ParticleList(): list<particle_t>()
ParticleList::ParticleList(): SmartList<Particle>()
{
}

ParticleList::~ParticleList()
{
}



point_t ParticleList::centerOfMassVelocity()
{
    point_t v = {{ 0, 0, 0 }};

    /*    for (particle_p i = begin(); i != end(); i++)
          v += i->v;*/

    SL_FOR_EACH
      (Particle,
       (*this),
       v += __iSLFE->v;
       );

    v/=size();

    return v;
}


void ParticleList::scaleVels(double temperature)
{
  bool many = size() > 1;
  
//   MSG_DEBUG("ParticleList::scaleVels", "size = " << size());  
  point_t v = centerOfMassVelocity();
//   MSG_DEBUG("ParticleList::scaleVels", "vCM = " << v);
  double vv = 0;
  double scale;

  
  /* No center of mass motion (COMM)... */
  SL_FOR_EACH
    (Particle,
     (*this),
//      MSG_DEBUG("ParticleList::scaleVels", "looping over P with v = " << i->v);
     // we can only remove COMM if there's more than one particle
     if(many) __iSLFE->v -= v;
     vv += __iSLFE->v.absSquare();
    );
//   MSG_DEBUG("ParticleList::scaleVels", "vv = " << vv);
  vv /= size();
//   MSG_DEBUG("ParticleList::scaleVels", "vv/size = " << vv);

  /* ... but the correct temperature. */
  scale = sqrt(SPACE_DIMS*temperature/vv);  
//   MSG_DEBUG("ParticleList::scaleVels", "scale = " << scale);

  /*    for (particle_p i = begin(); i != end(); i++)
        i->v *= scale;*/

  SL_FOR_EACH
    (Particle,
     (*this),
     __iSLFE->v *= scale;
    );
}
