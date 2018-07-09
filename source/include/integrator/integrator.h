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



#ifndef __INTEGRATOR_H
#define __INTEGRATOR_H

#include "node.h"
#include "smart_enum.h"
#include "particle.h"

#ifdef _OPENMP
  #include "omp.h"
#endif


using namespace std;

class Phase;
class Controller;

//---- Integrator ----

/*!
 * This is the base class for all integrators for an arbitraty degree of freedom.
 * Note that intgrators are responsible for registering species and the corresponding
 * degrees of freedom.
 */
class Integrator : public Node {
protected:
	/*!
	 * Species this integrator works for.
	 */
	string m_species;

#ifdef _OPENMP
  /*!
   * The offset to the copy-vectors in tag in \a Particle
   */
  vector<int> m_vec_offset;

  /*!
   * The offset to the proper place in the copy-vector for this \a Integrator 's copies
   */
  int m_vec_pos;

  /*!
   * Does this Integrator have to merge particle copies?
   */
  bool m_merge;
#endif

  /*!
   * Initialize the property list
   */
  void init();

	/*!
	 * Color this integrator works for.
	 */
	size_t m_colour;

	/*!
	 * Find the species with name \a species and return its color. Adds the species
	 * as new if it is not found. This simply calls the corresponding member
	 * of ManagerCell, which is protected and thus cannot be accessed by
	 * children of Integrator
	 * @param species Name of the species
	 */
	size_t getColourAndAdd(string species);

public:
	/*!
	 * Constructor
	 * @param controller Pointer to the controller this integrator belongs to
	 */
	Integrator(Controller *controller);

	/*!
	 * Destructor
	 */
	virtual ~Integrator();

	/*!
	 * Add species and set the correct name of this integrator.
	 */
	virtual void setup();

	/*!
	 * Called right before the simulation will start
	 */
	virtual void isAboutToStart() {
	}

	/*!
	 * If used forces are saved in a tag, protect and unprotect them as needed
	 */
	virtual void unprotect(size_t index) {
	}

	/*!
	 * Step 1 is called BEFORE update of the forces.
	 * Thus, advancing of positions should happen in step 1.
	 */
	virtual void integrateStep1() = 0;

#ifdef _OPENMP
  /*!
   * Returns the number of doubles the \a Integrator saves in a particle tag.
   */
  virtual int numCopyDoubles() = 0;

  /*!
   * Returns the offset to the copy tags of a \a Particle
   */
  virtual vector<int> &offsetToVec() {
    return m_vec_offset;
  }

  /*!
   * Returns the offset to the proper place in the copy-vector for this \a Integrator 's copies
   */
  virtual int &posInVec() {
    return m_vec_pos;
  }

  /*!
   * Degree of freedom for this \a Integrator
   */
  virtual string dofIntegr() = 0;

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index) = 0;

  /*!
   * Does this Integrator need to merge copies in Particle?
   */
  virtual bool &merge() {
    return m_merge;
  }
#endif

  /*!
   * Return the color this integrator works for.
   */
  size_t colour() const {
    return m_colour;
  }

	/*!
	 * Step 2 is called AFTER update of the forces.
	 */
	virtual void integrateStep2() = 0;

	/*!
	 * Calculate derived quantities, like, e.g., the temperature which can as well
	 * be calculated using the internal energy.
	 */
	virtual void deriveQuantities() {
	}

};

//---- Factories ----

class Integrator_Factory : public SmartEnum<Integrator_Factory> {
public:
	virtual Integrator *instantiate(Controller *controller) const = 0;

protected:
	Integrator_Factory(const string &name) :
		SmartEnum<Integrator_Factory>(name) {
	}
};

template <class T> class Integrator_Register : public Integrator_Factory {
public:
	Integrator_Register(const string &name) :
		Integrator_Factory(name) {
	}

	virtual Integrator *instantiate(Controller *controller) const;
};

//---- Inline functions ----

template <class T> inline Integrator *Integrator_Register<T>::instantiate(
		Controller *controller) const {
	return new T(controller);
}

#endif
