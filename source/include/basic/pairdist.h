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




#ifndef __PAIRDIST_H
#define __PAIRDIST_H

#include <vector>

#include <math.h>

//#include "misc.h"
#include "general.h"
#include "particle.h"
#include "data_format.h"
#include "smart_pointer.h"
#include "val_calculator.h"
#include "smart_list.h"

using namespace std;


//--- typedefs ----

/*!
 * Raw distance information
 */
struct dist_t {
  /*!
   * Distance in cartesian coordinates
   */
  point_t cartesian;

  /*!
   * Absolute value of the distance
   */
  double abs;

  /*!
   * Distance squared
   */
  double abs_square;

  /*!
   * Calculate \a abs and \a abs_square members
   */
  inline void calcAbs() {
    abs_square = 0;
    for (int i = 0; i < SPACE_DIMS; i++)
      abs_square += cartesian[i]*cartesian[i];
    abs = sqrt(abs_square);
  }
};



//---- Pairdist -----------------------------------------------------------

class ColourPair;

/*!
 * Information about the distance between two particle pairs
 */
class Pairdist
{
protected:
  /*!
   * The \a ColourPair for the two \a Particle s. Note: this is redundant
   * information.
   */
  ColourPair* m_cp;

// FIXME!!! m_distance is currently public because of
// offsetof(class Pairdist, m_distance) calls

//   dist_t m_distance;

  /*!
   * This \a Pairdist stores the distance between these two \a Particle 's
   */
  pair<Particle*, Particle*> m_particles;

  /*!
   * Does the force act on the first and/or the second particle?
   */
  pair<bool, bool> m_acts_on;

  /*!
   * Run the calculators for this \a Pairdist (stage 0)
   */
#ifndef _OPENMP
  virtual void initVals();
#endif

public:
  /*!
   * Default constructor. Only needed for vector<Pairdist>
   */
  Pairdist();

  /*!
   * Constructor
   * @param cp \a ColourPair for this \a Pairdist
   */
  Pairdist(ColourPair* cp);

  /*!
   * Copy constructor
   * @param pair \a Pairdist to copy
   */
  Pairdist(const Pairdist &pair);

  /*!
   * Constructor
   * @param cp \a ColourPair of this pair
   * @param dist Raw distance information (cartesian, abs)
   * @param first Pointer to the first particle in this pair
   * @param second Pointer to the second
   * @param ao_first Does the force act on the first particle?
   * @param ao_second Does the force act on the second particle?
   */
  Pairdist(ColourPair* cp, dist_t dist, Particle* first, Particle* second,
           bool ao_first = true, bool ao_second = true);

  /*!
   * Destructor
   */
  virtual ~Pairdist();

  /*!
   * Copy operator
   * @param pair \a Pairdist to copy
   */
  Pairdist &operator=(const Pairdist &pair);

  /*!
   * Set the properties of this \a Pairdist. Note that the \a ColourPair
   * has to be set beforhand, because \a set calls \a initVals which
   * calls the calculators. The information about the calculators, on the other
   * hand, is stored within the \a ColourPair.
   * @param v Raw distance information (cartesian, abs)
   * @param first Pointer to the first particle in this pair
   * @param second Pointer to the second
   * @param ao_first Does the force act on the first particle?
   * @param ao_second Does the force act on the second particle?
   */
  void set
    (dist_t v, Particle* first, Particle* second,
     bool ao_first = true, bool ao_second = true);

  /*!
   * Set the distance for new created pairs. Used by te \a VLYaoCreator
   * @param v Raw distance information (cartesian, abs)
   */
  void set(dist_t v);

  /*!
   * Set the \a ColourPair for this \a Pairdist
   * @param cp The \a ColourPair
   */
  void setCP(ColourPair* cp);

  /*!
   *
   */
  virtual ColourPair *cp() {
    return m_cp;
  }

  // CONVENTION 2: the following two functions will always return 'false' for a frozen particle
  // see also cell.cpp:CellLink::createDistances(..) for that

  /*!
   * Does this force act on the first particle?
   */
  inline bool actsOnFirst() const {
    return m_acts_on.first;
  }

  /*!
   * Does this force act on the second particle?
   */
  inline bool actsOnSecond() const {
    return m_acts_on.second;
  }

  /*!
   * Return the pointer to the first particle
   */
  inline Particle* firstPart()
  {
    return m_particles.first;
  }

  /*!
   * Return the pointer to the second particle
   */
  inline Particle* secondPart()
  {
    return m_particles.second;
  }

  /*!
   * Return the cartesian distance between the two particles
   */
  inline const point_t &cartesian() const {
    return m_distance.cartesian;
  }

  /*!
   * Return the absolute distance between the two particles
   */
  inline double abs() const {
    return m_distance.abs;
  }

  /*!
   * Return the square of the distance between the two particles
   */
  inline double absSquare() const {
    return m_distance.abs_square;
  }

  /*!
   * Run all standard calculators corresponding to stage \a stage for this \a Pairdist
   * @param stage Run calculators for this stage
   */
#ifndef _OPENMP
  void runCalculatorsForStage(size_t stage);
#else
  void runCalculatorsForStage(size_t stage, int thread_no);
#endif

  /*!
   * Run all standard calculators corresponding to stage \a stage for this \a Pairdist
   * @param stage Run calculators for this stage
   */
#ifndef _OPENMP
  void runBondedPairCalculators(size_t stage, size_t listIndex);
#else
  void runBondedPairCalculators(size_t stage, size_t listIndex, size_t thread_no);
#endif

  /*!
   * Run all additional calculators corresponding to stage \a stage for this \a Pairdist
   * @param stage Run calculators for this stage
   */
#ifndef _OPENMP
  void runBondedPairCalculators_0(size_t stage, size_t listIndex);
#else
  void runBondedPairCalculators_0(size_t stage, size_t listIndex, size_t thread_no);
#endif

  /*!
   * Run all additional calculators corresponding to stage \a stage for this \a Pairdist
   * @param stage Run calculators for this stage
   */
#ifndef _OPENMP
  void runCalculatorsForStage_0(size_t stage);
#else
  void runCalculatorsForStage_0(size_t stage, int thread_no);
#endif

  /*!
   * Update the distances of the particles including the absolute value
   */
  void calculateDistance();

  /*!
   * Update the cartesian distances of the the particles
   */
  void calculateCartDistance();

  /*!
   * Return the difference in cartesian coordinates
   * @param i x, y or z
   */
  inline double operator[](int i) const {
    assert(i >= 0 && i < SPACE_DIMS);
    if (i == SPACE_DIMS)
      return m_distance.abs_square;
    else
      return m_distance.cartesian[i];
  }

  /*!
   * The distance information
   */
  dist_t m_distance;

  /*!
   * Additional information to be stored in this \a Pairdist,
   * i.e., the reciprocal of the absolute value of the distance.
   */
  Data tag;

  SMARTLISTENTRY(Pairdist)
};

#endif
