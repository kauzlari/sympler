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



#include "pairdist.h"
#include "manager_cell.h"
#include "colour_pair.h"

using namespace std;


//---- Constructors/Destructor ----

// this constructor is currently needed for resizing Pairlists (01/13/05)
// immediately after the call, Pairdist::set(cp, ...) MUST be called (01/13/05)
Pairdist::Pairdist(): m_cp(NULL), m_particles(pair<Particle*, Particle*>(NULL, NULL)), tag()
{
//   throw gError("Pairdist::Pairdist(): don't call me!");
}

Pairdist::Pairdist(ColourPair* cp)
  : m_cp(cp), m_particles(pair<Particle*, Particle*>(NULL, NULL)), tag()
{
//   MSG_DEBUG
//     ("Pairdist::Pairdist(ColourPair*)",
//      "m_cp = (" << m_cp->firstSpecies() << ", " << m_cp->secondSpecies() << ")");

  tag.setFormatAndAlloc(&(m_cp->tagFormat()));
}


Pairdist::Pairdist(ColourPair* cp, dist_t v, Particle* first, Particle* second, bool ao_first, bool ao_second)
  : m_cp(cp), tag()
{
//   MSG_DEBUG
//     ("Pairdist::Pairdist(ColourPair*, ...)",
//      "m_cp = (" << m_cp->firstSpecies() << ", " << m_cp->secondSpecies() << ")");

  tag.setFormatAndAlloc(&(m_cp->tagFormat()));

  set(v, first, second, ao_first, ao_second);
}


Pairdist::Pairdist(const Pairdist &pair)
  : tag(pair.tag)
{
    m_cp = pair.m_cp;

//     MSG_DEBUG
//       ("Pairdist::Pairdist(Pairdist)",
//        "m_cp = (" << m_cp->firstSpecies() << ", " << m_cp->secondSpecies() << ")");

//     MSG_DEBUG
//       ("Pairdist::Pairdist(Pairdist)",
//        "tag.format()->rows() = " << tag.format()->rows() <<
//        ", tag.format().isNull() == " << tag.isNull());

    if (pair.m_particles.first) {
      m_distance = pair.m_distance;

      m_particles.first = pair.m_particles.first;
      m_particles.second = pair.m_particles.second;

      m_acts_on.first = pair.m_acts_on.first;
      m_acts_on.second = pair.m_acts_on.second;

//      tag.clear();
    }
/*      set(
	  pair.m_distance,
	  pair.m_particles.first,
	  pair.m_particles.second,
	  pair.m_acts_on.first,
	  pair.m_acts_on.second);*/
}


Pairdist::~Pairdist()
{
}

//---- Methods ----


Pairdist &Pairdist::operator=(const Pairdist &pair)
{
  tag = pair.tag;
  m_cp = pair.m_cp;

  if (pair.m_particles.first)
    set
      (pair.m_distance,
       pair.m_particles.first,
       pair.m_particles.second,
       pair.m_acts_on.first,
       pair.m_acts_on.second);
  return *this;
}

void Pairdist::setCP(ColourPair* cp)
{
  m_cp = cp;
  tag.setFormatAndAlloc(&(m_cp->tagFormat()));

}

void Pairdist::set(dist_t v, Particle* first, Particle* second, bool ao_first, bool ao_second)
{
  m_distance = v;

  m_particles.first = first;
  m_particles.second = second;

  m_acts_on.first = ao_first;
  m_acts_on.second = ao_second;
//   MSG_DEBUG("Pairdist::set", "calling initVals, CP = " << first->c << second->c);

//  initVals();
}


void Pairdist::set(dist_t v)
{/*
MSG_DEBUG("Pairdist::set", "m_distance BEFORE = " << m_distance.abs);*/
  m_distance = v;
// MSG_DEBUG("Pairdist::set", "m_distance AFTER = " << m_distance.abs);
}


#ifndef _OPENMP
void Pairdist::runCalculatorsForStage(size_t stage)
{
   FOR_EACH
      (vector<ValCalculator*>,
       m_cp->valCalculators(stage),
       (*__iFE)->compute(this);
      );
}

#else
void Pairdist::runCalculatorsForStage(size_t stage, int thread_no)
{
    FOR_EACH
      (vector<ValCalculator*>,
       m_cp->valCalculators(stage),
       // FIXME!: Are these checks really necessary?
       if (!(m_cp->valCalculators(stage).empty()))
	 (*__iFE)->compute(this, thread_no);
       );
    FOR_EACH
      (vector<ValCalculator*>,
       m_cp->valCalculatorParts(stage),
       if (!(m_cp->valCalculatorParts(stage).empty()))
	 (*__iFE)->compute(this, thread_no);
       );
}
#endif

#ifndef _OPENMP
void Pairdist::runBondedPairCalculators(size_t stage, size_t listIndex)
{
  FOR_EACH
      (vector<ValCalculator*>,
       m_cp->bondedValCalculators(stage, listIndex),
       (*__iFE)->compute(this);
      );
}

#else
void Pairdist::runBondedPairCalculators(size_t stage, size_t listIndex, size_t thread_no)
{
  vector<ValCalculator*>& vCs = m_cp->bondedValCalculators(stage, listIndex);
  FOR_EACH
    (vector<ValCalculator*>,
     vCs,
     if (!(vCs.empty()))
       (*__iFE)->compute(this, thread_no);
     );
}
#endif


#ifndef _OPENMP
void Pairdist::runBondedPairCalculators_0(size_t stage, size_t listIndex)
{
  FOR_EACH
      (vector<ValCalculator*>,
       m_cp->bondedValCalculators_0(stage, listIndex),
       (*__iFE)->compute(this);
      );
}

#else
void Pairdist::runBondedPairCalculators_0(size_t stage, size_t listIndex, size_t thread_no)
{
  vector<ValCalculator*>& vCs = m_cp->bondedValCalculators_0(stage, listIndex);
  FOR_EACH
    (vector<ValCalculator*>,
     vCs,
     if (!(vCs.empty()))
       (*__iFE)->compute(this, thread_no);
     );
}
#endif



#ifndef _OPENMP
void Pairdist::runCalculatorsForStage_0(size_t stage)
{
  FOR_EACH
      (vector<ValCalculator*>,
       m_cp->valCalculators_0(stage),
       (*__iFE)->compute(this);
      );
}
#else
void Pairdist::runCalculatorsForStage_0(size_t stage, int thread_no)
{
    FOR_EACH
        (vector<ValCalculator*>,
        m_cp->valCalculators_0(stage),
        (*__iFE)->compute(this, thread_no);
        );
    FOR_EACH
        (vector<ValCalculator*>,
        m_cp->valCalculatorParts_0(stage),
        (*__iFE)->compute(this, thread_no);
        );
}
#endif

void Pairdist::calculateDistance(void)
{
  m_distance.abs_square = 0;
  for (int _i = 0; _i < SPACE_DIMS; _i++)
    {
      m_distance.cartesian[_i] = m_particles.first->r[_i] - m_particles.second-> r[_i];
      m_distance.abs_square += m_distance.cartesian[_i]*m_distance.cartesian[_i];
    }
  /* Take care: The order of *j, *i defines the direction d	\
     is pointing to. */
  m_distance.abs = sqrt(m_distance.abs_square);
}

void Pairdist::calculateCartDistance(void)
{
  for (int _i = 0; _i < SPACE_DIMS; _i++)
    {
      m_distance.cartesian[_i] = m_particles.first->r[_i] - m_particles.second-> r[_i];
    }
  /* Take care: The order of *j, *i defines the direction d	\
     is pointing to. */
}

#ifndef _OPENMP
void Pairdist::initVals()
{
  runCalculatorsForStage(0);
}
#endif


