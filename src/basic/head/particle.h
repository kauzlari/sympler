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



#ifndef __PARTICLE_H
#define __PARTICLE_H

#include "smart_list.h"
#include "data_format.h"
#include "geometric_primitives.h"
// #include "particle_cache.h"

#ifdef ENABLE_PTHREADS
extern pthread_mutexattr_t g_mutex_attr;

#include "pthread.h"
#endif

// #define PCA_MAX_STAGE 2

/*--- Particle ---*/

#define FORCE_HIST_SIZE 2

class ParticleCache;
// class Pairdist;

/*!
 * Represent all the information that needs to be stored for each particle.
 * Furthermore, the \a tag field can used to expand this struture during run-time.
 */
class Particle
{
protected:
#ifdef ENABLE_PTHREADS
  pthread_mutex_t m_mutex; // Mutex for concurrent access
#endif

  /*!
   * Copy the information from particle \a p. Used by the copy
   * constructor and the copy operator.
   * @param p Particle to copy from
   */
  void copyFrom(const Particle &p) {
    /* Copy everything *except* mySlot, prev and next. */

    m_mass = p.m_mass;
    r = p.r;
    v = p.v;
    g = p.g;

/* TESTING whether needed*/
/*     displacement = p.displacement; */

    dt = p.dt;
    setColour(p.c);
    
    isFrozen = p.isFrozen;

    for (int i = 0; i < FORCE_HIST_SIZE; i++)
      force[i] = p.force[i];

    tag = p.tag;
  }

public:
  /*!
   * The particle mass. Currently (2011-12-06) not used and alway 1. Some 
   * Integrators also allow to set the mass. If you extend the 
   * usage of the mass here, take care for redundancy!!!
   */
  double m_mass;

  /*!
   * Stores the current position of the particle
   */
  point_t r;

  /*!
   * Stores the current velocity of the particle
   */
  point_t v;

  /*!
   * Stores the group the particle belongs to. The groups is
   * always the group of the cell the particle currently belongs to.
   */
  size_t g;

  /*!
   * The color, i.e., flavor of this particle
   */
  size_t c;                        /* colour */
  
  /*!
   * An attribute which states if the particle is free or frozen.
   */
    bool isFrozen;                 /* isFrozen */

  /*!
   * The displacement of the particle from its initial position. Used
   * for the calculation of the self-diffusion coefficient.
   */
/* TESTING whether needed*/
/*   point_t displacement; */

  /*!
   * Used for collision detection. so the integrator knows which dt to use.
   * This means, whenever the particle hits a wall the value of \a dt is modified to
   * the remaining time.
   */
  double dt; 

  /*!
   * The current force on the particle plus the last force from the last timestep.
   * Needed for the Velocity-Verlet algorithm.
   */
  point_t force[FORCE_HIST_SIZE];

  /*!
   * Additional data, like local densities, temperature, etc.
   */
  Data tag;

  /*!
   * A list of neighbouring particles
   */
//   vector< list<size_t> > m_pairs;

  /*!
   * Constructor
   */ 
  Particle(): m_mass(1.), c((size_t)HUGE_VAL) {
#ifdef ENABLE_PTHREADS
    pthread_mutex_init(&m_mutex, &g_mutex_attr);
#endif
//     displacement.assign(0);
    isFrozen = 0;
    clear();
  }

  /*!
   * Constructor
   * @param colour The color of the newly created particle
   */
    Particle(size_t colour): m_mass(1.), c((size_t)HUGE_VAL) {
#ifdef ENABLE_PTHREADS
    pthread_mutex_init(&m_mutex, &g_mutex_attr);
#endif
    setColour(colour);
//     displacement.assign(0);
    isFrozen = 0;
    clear();
  }

  /*!
   * Copy constructor
   * @param p The particle to copy from
   */
  Particle(const Particle& p) {
#ifdef ENABLE_PTHREADS
    pthread_mutex_init(&m_mutex, &g_mutex_attr);
#endif
    copyFrom(p);
  } 
       
  /*!
   * Destructor
   */
  ~Particle() {
#ifdef ENABLE_PTHREADS
    pthread_mutex_destroy(&m_mutex);
#endif
  }


  /*!
   * Assignment operator
   * @param p The particle to copy from
   */
  Particle &operator=(const Particle &p) {
    copyFrom(p);

    return *this;
  }

  /*! 
   * Set the particles' colour
   * @param newC The new color
   */
  void setColour(size_t newC) {
    if (tag.isNull() || c != newC) {
      if(!tag.isNull())
        tag.release();

      c = newC;
//      assert(c < s_tag_format.size());

      /* If this is not the case, c is HUGE_VAL which means
         not assigned. */
      if (c < s_tag_format.size()) {
          tag = Data(s_tag_format[c]);
      }
    }
  }

  /*!
   * Set the force, stress tensor, ... to zero
   * @param force_index The force index to clear
   */
  void clear(size_t force_index) {
    force[force_index].assign(0);
//    stress.assign(0); 
/* TESTING WHETHER NEEDED */
/*     displacement.assign(0); */
    //    virial_sum = 0;
  }

  void clear() {
    for(size_t i = 0; i < FORCE_HIST_SIZE; ++i)
      force[i].assign(0);
/* TESTING WHETHER NEEDED */
/*     displacement.assign(0); */
  }

#ifdef ENABLE_PTHREADS
  void lock() {
    pthread_mutex_lock(&m_mutex);
  }
  void unlock() {
    pthread_mutex_unlock(&m_mutex);
  }
#endif

  /*!
 * Properties that are recalculated for every timestep 
 * from information given only within the local particle,
 * i.e., a cache.
 * outermost vector is for colours and holds vectors for the stages, which 
 * hold the ParticleCache vectors 
   */
  static vector<vector<vector<ParticleCache*> > > s_cached_properties;
  

  /*!
   * This is a helper for building \a s_cached_properties . First, the \a ParticleCahce s
   * are registered here and afterwards sorted by stages into \a s_cached_properties .
   */
  static vector<ParticleCache*> s_cached_flat_properties;
  
  /*!
   * The biggest stage occuring in \a s_cached_properties . This is determined 
   * during runtime in Particle::sortStages()
   */
  static size_t s_maxStage;

  /*!
   * Same as above, but the properties are computed at another stage in the timestep, 
   * determined by the \a Controller 
   */
  static vector<vector<vector<ParticleCache*> > > s_cached_properties_0;
  

  /*!
   * This is a helper for building \a s_cached_properties_0 . First, the \a ParticleCache s
   * are registered here and afterwards sorted by stages into \a s_cached_properties_0 .
   */
  static vector<ParticleCache*> s_cached_flat_properties_0;
  
  /*!
   * The biggest stage occuring in \a s_cached_properties_0 . This is determined 
   * during runtime in Particle::sortStages_0()
   */
  static size_t s_maxStage_0;

  /*!
   * Description of the field \a tag.
   */
  static vector<Data> s_tag_format;
  
  /*!
   * Return the format description for a certain color
   * @param colour The color
   */
  static Data &tagFormat(size_t colour) {
    return s_tag_format[colour];
  }

  /*!
   * Register a new cache function in \a s_cached_flat_properties.
   * Will check if the cache function already exists,
   * if it does return the existing cache function.
   * @param c The new cache function
   */
  static ParticleCache* registerCache(ParticleCache *c);

  /*!
   * Register a new cache function in \a s_cached_flat_properties_0.
   * Will check if the cache function already exists,
   * if it does, return the existing cache function.
   * @param c The new cache function
   */
  static ParticleCache* registerCache_0(ParticleCache *c);

    /*!
   * Sort the \a ParticleCache s from \a s_cached_flat_properties into 
   * \a s_cached_properties according to their stages
     */
  static void sortStages();

      /*!
   * Sort the \a ParticleCache s from \a s_cached_flat_properties_0 into 
   * \a s_cached_properties_0 according to their stages
       */
  static void sortStages_0();

  SMARTLISTENTRY(Particle)
};

#endif
