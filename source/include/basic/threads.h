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


#ifndef __THREADS_H
#define __THREADS_H

#include "colour_pair.h"

/*! 
 * The code her is NOT parallel. We had some pthread versions here,
 * but they are removed now because they were inefficient, and we 
 * would nowadays use OpenMP or similar
 */

#define LL_FOR_EACH__PARALLEL(type, first, size, add_data, code)  \
{                                                                 \
  /*const size_t thread_no = 0;*/					\
  void *data = add_data;						      \
  type *next, *i;							\
  i = first;                                                      \
  while (i) {                                                     \
    next = i->next;                                               \
    code                                                          \
    i = next;                                                     \
  }                                                               \
} while(0)


/* #define FOR_EACH__PARALLEL(type, stlcont, add_data, code)	\ */
#define FOR_EACH__PARALLEL(type, stlcont, code)     \
{                                                       \
  /*const size_t thread_no = 0;*/			\
  /* void *data = add_data;*/				\
  FOR_EACH(type, stlcont, code);                        \
} while(0)


#define FOR_EACH_PAIR__PARALLEL(parent_type, colourpair, code)          \
{                                                                       \
  parent_type *self;							\
  self = this;								\
  FOR_EACH_PAIR(colourpair, code);                                      \
} while(0)


/* --- parallel version of FOR_EACH_PARTICLE_* --- */
#define FOR_EACH_FREE_PARTICLE_C__PARALLEL(phase, col, add_data, code)  \
{                                                                       \
  void *data = add_data;                                                \
  if (col == ALL_COLOURS) {                                             \
    size_t __ncols = phase->nColours();					\
    for(size_t __c = 0; __c < __ncols; ++__c) {                     \
      LL_FOR_EACH__PARALLEL                                             \
        (Particle,                                                      \
         phase->particles(__c).first(),                                   \
         phase->particles(__c).size(),                                    \
         add_data,                                                      \
         code);                                                         \
    }                                                                   \
  } else {                                                              \
      LL_FOR_EACH__PARALLEL                                             \
        (Particle,                                                      \
         phase->particles(col).first(),                                   \
         phase->particles(col).size(),                                    \
         add_data,                                                      \
         code);                                                         \
  }                                                                     \
} while(0)


#define FOR_EACH_FREE_PARTICLE__PARALLEL(phase, col, add_data, code)  \
{                                                                     \
  void *data = add_data;                                              \
  for(size_t c = 0; c < phase->nColours(); ++c) {                     \
    LL_FOR_EACH__PARALLEL                                             \
      (Particle,                                                      \
       phase->particles(c).first(),                                   \
       phase->particles(c).size(),                                    \
       add_data,                                                      \
       code);                                                         \
  }                                                                   \
} while(0)



#define FOR_EACH_PARTICLE_C__PARALLEL(phase, col, add_data, code)       \
{                                                                       \
  void *data = add_data;                                                \
  if (col == ALL_COLOURS) {                                             \
    for(size_t c = 0; c < phase->nColours(); ++c) {                     \
      LL_FOR_EACH__PARALLEL                                             \
        (Particle,                                                      \
         phase->particles(c).first(),                                   \
         phase->particles(c).size(),                                    \
         add_data,                                                      \
         code);                                                         \
      LL_FOR_EACH__PARALLEL                                             \
        (Particle,                                                      \
         phase->frozenParticles(c).first(),                             \
         phase->frozenParticles(c).size(),                              \
         add_data,                                                      \
         code);                                                         \
    }                                                                   \
  } else {                                                              \
    size_t c = col;                                                     \
      LL_FOR_EACH__PARALLEL                                             \
        (Particle,                                                      \
         phase->particles(c).first(),                                   \
         phase->particles(c).size(),                                    \
         add_data,                                                      \
         code);                                                         \
      LL_FOR_EACH__PARALLEL                                             \
        (Particle,                                                      \
         phase->frozenParticles(c).first(),                             \
         phase->frozenParticles(c).size(),                              \
         add_data,                                                      \
         code);                                                         \
  }                                                                     \
} while(0)


#define FOR_EACH_PARTICLE__PARALLEL(phase, col, add_data, code) \
{                                                               \
  void *data = add_data;                                        \
  for(size_t c = 0; c < phase->nColours(); ++c) {               \
    LL_FOR_EACH__PARALLEL                                       \
      (Particle,                                                \
       phase->particles(c).first(),                             \
       phase->particles(c).size(),                              \
       add_data,                                                \
       code);                                                   \
    LL_FOR_EACH__PARALLEL                                       \
      (Particle,                                                \
       phase->frozenParticles(c).first(),                       \
       phase->frozenParticles(c).size(),                        \
       add_data,                                                \
       code);                                                   \
  }                                                             \
} while(0)


#endif
