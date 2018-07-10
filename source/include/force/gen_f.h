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



#ifndef __GEN_F_H
#define __GEN_F_H

#include "node.h"
#include "phase.h"
#include "smart_enum.h"
#include "pairdist.h"

#ifdef _OPENMP
  #include "omp.h"
#endif

using namespace std;

/*---- GenF ----*/

struct rf_t {
  double xy;
  double yz;
  double zx;
};


class Simulation;
class ColourPair;
#ifdef _OPENMP
class Integrator;
#endif

/*!
 * Base class for all forces
 */
class GenF: public Node
{
public:
  /*!
   * Default constructor
   */
  GenF();

  /*!
   * Constructor that should be used
   */
  GenF(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~GenF();

  /*!
   * Compute all forces acting between pairs and store them in the \a force_index history slot
   * @param force_index The history slot to use
   */
#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index) = 0;
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no) = 0;

  /*!
   * Merges the copies from the copy-tags in \a Particle ad adds them to the real tag.
   */
//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) = 0;
#endif

  /*!
   * Compute all forces acting on a particle and store them in the \a force_index history slot
   * @param force_index The history slot to use
   */
  virtual void computeForces(Particle* part, int force_index) = 0;

  /*!
   * Compute all forces and store them in the \a force_index history slot
   * @param force_index The history slot to use
   */
  virtual void computeForces(int force_index) = 0;

  /*!
   * \a invalidate tells the force to erase all the temporary information,
   * i.e. force factors, etc.
   */
  inline virtual void invalidate() {
  }

   /*!
   * Return the name of this force function
   */
  const string &name() const {
    return m_force_name;
  }

  /*!
   * Inititialize the property list.
   */
  void init();

  /*!
   * Return the tag that tells the Simulation whether it is a particle or pair force, to be added in the proper forces-vector
   */
  virtual bool isParticleForce() {
    return m_is_particle_force;
  }

  /*!
   * Return the tag that tells the Simulation whether it is a particle or pair force, to be added in the proper forces-vector
   */
  virtual bool isPairForce() {
    return m_is_pair_force;
  }

#ifdef _OPENMP
  /*!
   * Returns the offset to the copy tags of a \a Particle
   */
  virtual vector<pair<int, int> > &vecOffset() {
    return m_offsetToVec;
  }

  /*!
   * Returns the offset to the proper place in the copy-vector for this \a Force 's copies
   */
  virtual pair<int, int> &vecPos() {
    return m_posInVec;
  }

  /*!
   * Sets the offset to the copy-slots in \a Particle 's tag
   */
  virtual void setForceSlots(Integrator* intr, int thread_no) = 0;

#endif

  virtual size_t fCol() {
    return m_colour;
  }

protected:

  /*!
   * The name of the weigting function.
   */
  string m_force_name;

#ifdef _OPENMP
  /*!
   * The offset to the copy-vectors in tag in \a Particle
   */
  vector<pair<int, int> > m_offsetToVec;

  /*!
   * The offset to the proper place in the copy-vector for this \a Force 's copies
   */
  pair<int, int> m_posInVec;
#endif

  /*!
   * Tag that tells the Simulation whether it is a particle force, to be added in the proper forces-vector
   */
  bool m_is_particle_force;

  /*!
   * Tag that tells the Simulation whether it is a pair force, to be added in the proper forces-vector
   */
  bool m_is_pair_force;

  /*!
   * The colour index that corresponds to the species \a m_species
   */
  size_t m_colour;
};


/*---- Factory ----*/

class GenFType : public SmartEnum<class GenFType>
{
 public:
  virtual class GenF* makeGenF(Simulation *simulation) const = 0;
 protected:
  GenFType(const string& name)
	: SmartEnum<GenFType>(name){}
};

template <class GenFConcr>
class GenFTypeConcr : public GenFType
{
 public:
  virtual class GenF* makeGenF(Simulation *simulation) const;
  GenFTypeConcr(const string& name)
	: GenFType(name) {}
};

//---- inline functions ---------------------------------------------------

template <class GenFConcr>
inline GenF* GenFTypeConcr<GenFConcr>::makeGenF(Simulation *simulation) const
{
  return new GenFConcr(simulation);
}



//---- Functions ----

group_t group_intersection(const group_t &gr1, const group_t &gr2);

#endif


