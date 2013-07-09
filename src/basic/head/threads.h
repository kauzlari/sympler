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


#ifndef __THREADS_H
#define __THREADS_H

#include "config.h"
#include "colour_pair.h"

#ifdef ENABLE_PTHREADS

#include <pthread.h>

/* --- Thread macros --- */

#define THREAD_CHUNK_GRAIN 10

/* ENABLE_PTHREADS case

The current implementation of the parallelization works as follows:

The method operator() of the __Handle__Thread class is the actual thread.
This is done this way so we can transparently define parallelization by
just redefining the macro used.

The loop is parallelized by splitting it up into THREAD_CHUNK_GRAIN parts.
These are then executed in parallel. For example for two threads

  eeeeeeeeeeee|eeeeeeeeeeee|eee...
    Thread 1     Thread 2     Thread 1 ...

Each thread is skipping the parts it is not responsible for by calling
iterator->next.

*/



/* Base class */

struct BasicThread {
  virtual void operator()(size_t thread_no) = 0;
};

#if 0
/* --- for each loop for linked lists --- */
#define LL_FOR_EACH__PARALLEL(type, first, size, add_data, code)        \
{                                                                       \
  struct __Handle__Thread: public BasicThread {                         \
    type *cur;                                                          \
    void *data;                                                         \
    pthread_mutex_t mutex;                                              \
    size_t chunk_size;                                                  \
                                                                        \
    __Handle__Thread(type *c, size_t cs, void *d)                       \
      : cur(c), chunk_size(cs), data(d) {                               \
      pthread_mutex_init(&mutex, &g_mutex_attr);                        \
    }                                                                   \
    virtual ~__Handle__Thread() {                                       \
      pthread_mutex_destroy(&mutex);                                    \
    }                                                                   \
    virtual void operator()(size_t thread_no) {                         \
      /* We need two because ->next might change                        \
         between calls. */                                              \
      type *next, *i;                                                \
                                                                        \
      size_t debug_n_calls = 0;                                         \
                                                                        \
      pthread_mutex_lock(&mutex);                                       \
      i = cur;                                                       \
      for (size_t __i = 0; __i < chunk_size && cur; ++__i)              \
        cur = cur->next;                                                \
      pthread_mutex_unlock(&mutex);                                     \
                                                                        \
      while (i) {                                                    \
        for (size_t __i = 0; __i < chunk_size && i; ++__i) {         \
          debug_n_calls++;                                              \
          next = i->next;                                            \
          code                                                          \
                                                                        \
            i = next;                                                \
        }                                                               \
                                                                        \
        pthread_mutex_lock(&mutex);                                     \
        i = cur;                                                     \
        for (size_t __i = 0; __i < chunk_size && cur; ++__i)            \
          cur = cur->next;                                              \
        pthread_mutex_unlock(&mutex);                                   \
      }                                                                 \
                                                                        \
      /*    pthread_mutex_lock(&mutex);*/                               \
      /*    MSG_DEBUG(string(#name) + "__Thread::operator()", "thread_no = " << thread_no << ", n_calls = " << debug_n_calls); */ \
      /*    pthread_mutex_unlock(&mutex);*/                             \
    }                                                                   \
  };                                                                    \
                                                                        \
  __Handle__Thread *t;                                                  \
  size_t chunk_size;                                                    \
                                                                        \
  chunk_size = size/(global::n_threads*THREAD_CHUNK_GRAIN);             \
  /*chunk_size = 1;*/                                                   \
  if (chunk_size < 1)                                                   \
    chunk_size = 1;                                                     \
  /*  MSG_DEBUG(string(#name) + "->RUN_FOR_EACH_THREAD_LL", "chunk_size = " << chunk_size);*/ \
                                                                        \
  ForkThreads                                                           \
    (t = new __Handle__Thread(first, chunk_size, add_data));            \
  delete t;                                                             \
} while(0)
#endif


/* --- for each loop for stl lists/vectors --- */
#define FOR_EACH__PARALLEL(type, stlcont, add_data, code)               \
{                                                                       \
  struct __Handle__Thread: public BasicThread {                         \
    type::iterator cur;                                                 \
    type::iterator end_i;                                               \
    void *data;                                                         \
    pthread_mutex_t mutex;                                              \
    size_t chunk_size;                                                  \
                                                                        \
    __Handle__Thread(type::iterator bi, type::iterator ei, size_t cs, void *d) \
      : cur(bi), end_i(ei), chunk_size(cs), data(d) {                   \
      pthread_mutex_init(&mutex, &g_mutex_attr);                        \
    }                                                                   \
    ~__Handle__Thread() {                                               \
      pthread_mutex_destroy(&mutex);                                    \
    }                                                                   \
    void operator()(size_t thread_no) {                                 \
      type::iterator i;                                                 \
                                                                        \
      size_t debug_n_calls = 0;                                         \
                                                                        \
      pthread_mutex_lock(&mutex);                                       \
      i = cur;                                                          \
      if (cur != end_i) {                                               \
        for (size_t __i = 0; __i < chunk_size && cur != end_i; ++__i)   \
          cur++;                                                        \
      }                                                                 \
      pthread_mutex_unlock(&mutex);                                     \
                                                                        \
      while (i != end_i) {                                              \
        for (size_t __i = 0; __i < chunk_size && i != end_i; ++__i) {   \
          debug_n_calls++;                                              \
          code                                                          \
                                                                        \
            i++;                                                        \
        }                                                               \
                                                                        \
        pthread_mutex_lock(&mutex);                                     \
                                                                        \
        i = cur;                                                        \
        if (cur != end_i) {                                             \
          for (size_t __i = 0; __i < chunk_size && cur != end_i; ++__i) \
            cur++;                                                      \
        }                                                               \
        pthread_mutex_unlock(&mutex);                                   \
      }                                                                 \
                                                                        \
      /*    pthread_mutex_lock(&mutex);*/                               \
      /*    MSG_DEBUG(string(#name) + "__Thread::operator()", "thread_no = " << thread_no << ", n_calls = " << debug_n_calls); */ \
      /*    pthread_mutex_unlock(&mutex);*/                             \
    }                                                                   \
  };                                                                    \
                                                                        \
  __Handle__Thread *t;                                                  \
  size_t chunk_size;                                                    \
                                                                        \
  chunk_size = stlcont.size()/(global::n_threads*THREAD_CHUNK_GRAIN);   \
  if (chunk_size < 1)                                                   \
    chunk_size = 1;                                                     \
  /*  MSG_DEBUG(string(#name) + "->RUN_FOR_EACH_THREAD_STL", "chunk_size = " << chunk_size);*/ \
                                                                        \
  ForkThreads                                                           \
    (t = new __Handle__Thread(stlcont.begin(), stlcont.end(), chunk_size, add_data)); \
  delete t;                                                             \
} while(0)


/* --- loop over pairs */
/*FIXME!! currently(2006/05/12), we use the non-parallel version here, because none of the parallel solutions below works satisfactorily. But since the pair operations are the most expensive ones, the current parallel verion of the code (e.g. parallel loop over particles) gives almost no advantage*/
#define FOR_EACH_PAIR__PARALLEL(parent_type, colourpair, code)          \
{                                                                       \
  parent_type *self;							\
  self = this;								\
  FOR_EACH_PAIR(colourpair, code);                                      \
} while(0)



/* --- loop over pairs - BAD VERSION--- */
#if 0
#define FOR_EACH_PAIR__PARALLEL(parent_type, colourpair, code)          \
{\
  for(size_t __slot = 0; __slot < global::n_threads; ++__slot)	\
    LL_FOR_EACH__PARALLEL(Pairdist,				\
			  colourpair->pairs(__slot).first(),	\
			  colourpair->pairs(__slot).size(),\
			  NULL,\
			  Pairdist* pair = i;	\
			  code);		\
} while(0)
#endif


/* --- loop over pairs - THE BEST SOLUTION WE HAVE, BUT BAD TOO--- */
#if 0
#define FOR_EACH_PAIR__PARALLEL(parent_type, colourpair, code)          \
{                                                                       \
  struct __Handle__Thread: public BasicThread {                         \
    ColourPair *cp;                                                     \
    Pairdist* cur;							\
    size_t cur_slot;                                                    \
    parent_type *self;							\
    pthread_mutex_t mutex;                                              \
    size_t chunk_size;                                                  \
                                                                        \
    __Handle__Thread(ColourPair *c, size_t cs, parent_type *s)          \
      : cp(c), cur_slot(0), chunk_size(cs) , self(s)  {			\
      cur = cp->pairs(0).first();                                       \
      pthread_mutex_init(&mutex, &g_mutex_attr);                        \
    }                                                                   \
    ~__Handle__Thread() {                                               \
      pthread_mutex_destroy(&mutex);                                    \
    }                                                                   \
    void operator()(size_t thread_no) {                                 \
      Pairdist *pair, *last_pair;					\
                                                                        \
      size_t debug_n_calls = 0;                                         \
                                                                        \
      /**/ pthread_mutex_lock(&mutex);					\
      pair = cur;                                                       \
      if (cur != NULL) {						\
        for (size_t __i = 0; __i < chunk_size && cur != NULL; ++__i)	\
          cur = cur->next;						\
      }                                                                 \
      last_pair = cur;							\
      if (cur == NULL && cur_slot < 2*global::n_threads) {		\
        cur_slot++;                                                     \
        if (cur_slot < 2*global::n_threads) {                           \
          cur = cp->pairs(cur_slot).first();                            \
        }                                                               \
      }                                                                 \
      /**/ pthread_mutex_unlock(&mutex);				\
                                                                        \
      while (pair != NULL) {						\
        for (; pair != last_pair; pair = pair->next) {			\
          debug_n_calls++;                                              \
          code                                                          \
        }                                                               \
                                                                        \
        /**/ pthread_mutex_lock(&mutex);				\
        pair = cur;                                                     \
        if (cur != NULL) {						\
          for (size_t __i = 0; __i < chunk_size && cur != NULL; ++__i)	\
            cur = cur->next;						\
        }                                                               \
	last_pair = cur;						\
        if (cur == NULL && cur_slot < 2*global::n_threads) {		\
          cur_slot++;                                                   \
          if (cur_slot < 2*global::n_threads) {                         \
            cur = cp->pairs(cur_slot).first();                          \
          }                                                             \
        }                                                               \
        /**/ pthread_mutex_unlock(&mutex);				\
      }                                                                 \
                                                                        \
      /*    pthread_mutex_lock(&mutex);*/                               \
      /*    MSG_DEBUG(string(#name) + "__Thread::operator()", "thread_no = " << thread_no << ", n_calls = " << debug_n_calls); */ \
      /*    pthread_mutex_unlock(&mutex);*/                             \
    }                                                                   \
  };                                                                    \
                                                                        \
  __Handle__Thread *t;                                                  \
  size_t chunk_size = 0;                                                \
                                                                        \
  for (int i = 0; i < 2*global::n_threads; ++i) {                       \
    chunk_size = max(chunk_size, colourpair->pairs(i).size());          \
  }                                                                     \
                                                                        \
  chunk_size /= global::n_threads*THREAD_CHUNK_GRAIN;                   \
  if (chunk_size < 1)                                                   \
    chunk_size = 1;                                                     \
									\
  MSG_DEBUG								\
    (string(#parent_type) + "->FOR_EACH_PAIR__PARALLEL",		\
     "chunk_size = " << chunk_size					\
     );									\
                                                                        \
  ForkThreads                                                           \
    (t = new __Handle__Thread(colourpair, chunk_size, this));           \
  delete t;                                                             \
} while(0)
#endif



typedef void *(*thread_func_t)(void*);

struct thread_info_t {
  pthread_t id;
  size_t no;

  void *handler;
};


extern pthread_mutexattr_t g_mutex_attr;

void ForkThreads(BasicThread *handler);
void ForkThreadsWithFunction(thread_func_t function, void *handler);

// #else of #ifdef ENABLE_PTHREADS
#else
#define LL_FOR_EACH__PARALLEL(type, first, size, add_data, code)  \
{                                                                 \
  const size_t thread_no = 0;						  \
  void *data = add_data;						      \
  type *next, *i;							\
  i = first;                                                      \
  while (i) {                                                     \
    next = i->next;                                               \
    code                                                          \
    i = next;                                                     \
  }                                                               \
} while(0)


#define FOR_EACH__PARALLEL(type, stlcont, add_data, code)     \
{                                                             \
  const size_t thread_no = 0;					      \
  void *data = add_data;						      \
  FOR_EACH(type, stlcont, code);                              \
} while(0)


#define FOR_EACH_PAIR__PARALLEL(parent_type, colourpair, code)          \
{                                                                       \
  parent_type *self;							\
  self = this;								\
  FOR_EACH_PAIR(colourpair, code);                                      \
} while(0)

#endif


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
