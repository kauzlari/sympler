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


#ifndef __PARTICLE_CREATOR_FREE_PCALC_H
#define __PARTICLE_CREATOR_FREE_PCALC_H

#include "pc_free.h"

//---- Forward declarations ----


//---- Classes ----

class ParticleCreatorFree;
class ParticleCalculator;



struct dof_info_t {
  size_t offset;
  ParticleCalculator *pc;
//   string name;
};


class ParticleCalculator: public calculator
{
  protected:
    Particle *m_particle;
    size_t m_colour;
    point_t *m_box_size;

    virtual void AnalyzeId(const string &id);
	
  public:
    ParticleCalculator(): calculator(), m_particle(NULL), m_colour(0), m_box_size(NULL) {
      throw gError("ParticleCalculator::ParticleCalculator", "Do not use.");
    }
    ParticleCalculator(Particle *particle, size_t colour, point_t *box_size)
  : calculator(), m_particle(particle), m_colour(colour), m_box_size(box_size) {
  }
  ParticleCalculator(const ParticleCalculator& old)
  : calculator(), m_particle(old.m_particle), m_colour(old.m_colour),
  m_box_size(old.m_box_size) {
    FromString(old.text);
  }
  ParticleCalculator(Particle *particle, size_t colour, point_t *box_size, const string& old)
  : calculator(), m_particle(particle), m_colour(colour), m_box_size(box_size) {
    FromString(old);
  }
};


/*!
*\a ParticleCreatorFree using old-style \a ParticleCalculators for setting of initial values.
*/
class ParticleCreatorFreePCalc: public ParticleCreatorFree
{
  protected:
  
    /* Calculators. */
    vector<ParticleCalculator*> m_poss[SPACE_DIMS], m_vels[SPACE_DIMS], m_internal_energy;
    /* Position information used by calculators. */
    vector<Particle> __m_particle;

    string m_str_poss[SPACE_DIMS], m_str_vels[SPACE_DIMS];
    //   /* List of internal degrees of freedom to be set for each colour. */
    vector< vector<dof_info_t> > m_internal_dofs;

    void init();
    void initTransform();

    /* Transform applies the functions given by posX...Z and velX...Z */
    void transformPos(Particle &p) {
      if(__m_particle.size() == 1)
      {

        __m_particle[0] = p;

        for (int d = 0; d < SPACE_DIMS; d++){

          p.r[d] = m_poss[d][0]->eval();
}
  
        for (vector<dof_info_t>::iterator i = m_internal_dofs[0].begin();
             i != m_internal_dofs[0].end(); ++i) {
               p.tag.doubleByOffset(i->offset) = i->pc->eval();
             }
      }
      else
      {
        __m_particle[p.c] = p;
  
        for (int d = 0; d < SPACE_DIMS; d++)
          p.r[d] = m_poss[d][p.c]->eval();
  
        for (vector<dof_info_t>::iterator i = m_internal_dofs[p.c].begin();
             i != m_internal_dofs[p.c].end(); ++i) {
               p.tag.doubleByOffset(i->offset) = i->pc->eval();
        
             }
      }
    }

    void transformVel(Particle &p) {
      if(__m_particle.size() == 1)
      {
        __m_particle[0] = p;
  
        for (int d = 0; d < SPACE_DIMS; d++) {
          p.v[d] = m_vels[d][0]->eval();
        }
      }
      else
      {
        __m_particle[p.c] = p;
  
        for (int d = 0; d < SPACE_DIMS; d++) {
          p.v[d] = m_vels[d][p.c]->eval();
        }
      }
    }

	
  public:
    ParticleCreatorFreePCalc();
    ParticleCreatorFreePCalc(Boundary *boundary);
    virtual ~ParticleCreatorFreePCalc();

    virtual void setup();
/*    virtual void flushParticles();
    virtual void flushParticles(Particle** first_p);
    virtual ostream &write(ostream &s, int shift = 0);*/
};


/*
//---- Factories ----

class ParticleCreatorFree_Factory: public SmartEnum<ParticleCreatorFree_Factory>
{
  public:
    virtual ParticleCreatorFree *instantiate(Boundary *boundary) const = 0;

  protected:
    ParticleCreatorFree_Factory(const string &name)
  : SmartEnum<ParticleCreatorFree_Factory>(name) { }
};


template <class T>
    class ParticleCreatorFree_Register: public ParticleCreatorFree_Factory
{
  public:
    ParticleCreatorFree_Register(const string &name)
  : ParticleCreatorFree_Factory(name) { }
	
    virtual ParticleCreatorFree *instantiate(Boundary *boundary) const;
};



//---- Inline functions ----

template <class T>
    inline ParticleCreatorFree *ParticleCreatorFree_Register<T>::instantiate(Boundary *boundary) const
{
  return new T(boundary);
}
*/

#endif
