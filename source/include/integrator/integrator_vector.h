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


#ifndef __INTEGRATOR_VECTOR_H
#define __INTEGRATOR_VECTOR_H

#include <math.h>

// #include "function.h"
#include "particle.h"
#include "integrator.h"

using namespace std;

class GenF;
class Phase;
class Controller;

/*!
 * Simple Euler integration of a vector field
 */
class IntegratorVector: public Integrator
{
  protected:
  /*!
   * The name of the vector field
   */
    string m_vector_name;

  /*!
     * The symbol (short name) of the vector field
   */
    string m_vector_symbol;

  /*!
     * The tag offset of the vector field
   */
    size_t m_vector_offset;

  /*!
     * The tag offsets of the force on this vector field. Even though this integrator
     *needs only one, we let it use as much as the most sophisticated integrator. This is
     *currently (2006/01/11) FORCE_HIST_SIZE=2
   */
    size_t m_force_offset[FORCE_HIST_SIZE];

  /*!
   * The indices of the force attribute of the tag of the force on this vector field.    */
    size_t m_fAttr_index[FORCE_HIST_SIZE];

  /*!
     * Time step
   */
    double m_dt;

  /*!
     * Initialize the property list
   */
    void init();

  public:
  /*!
   * Constructor
   * @param controller Pointer to the \a Controller object this \a Integrator belongs to
   */
    IntegratorVector(Controller *controller);

  /*!
     * Destructor
   */
    virtual ~IntegratorVector();

  /*!
     * Register the field and the force of the field with the \a Particle
   */
    virtual void setup();

#ifdef _OPENMP
  /*!
   * Returns the number of doubles this Integrator saves in a particle tag.
   */
    virtual int numCopyDoubles() {
      return SPACE_DIMS;
    }

    virtual string dofIntegr();

  /*!
   * Merge the copies at the end of every timestep
   */
  virtual void mergeCopies(Particle* p, int thread_no, int force_index);

#endif

  /*!
     * Return the name of the vector field
   */
    virtual string vectorName() const {
      return m_vector_name;
    }

  /*!
     * If used forces are saved in a tag, protect and unprotect them as needed
   */
    virtual void unprotect(size_t index);

  /*!
     * Called right before the simulation will start
   */
    virtual void isAboutToStart();

  /*!
     * Return the tag offset of the vector field
   */
    virtual size_t vectorOffset() const {
      return m_vector_offset;
    }

  /*!
     * Integrate the field
   */
    virtual void integrateStep1();

  /*!
     * Does nothing
   */
    virtual void integrateStep2();
};

#endif
