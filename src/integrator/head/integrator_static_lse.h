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



#ifndef INTEGRATOR_STATIC_LSE_H_
#define INTEGRATOR_STATIC_LSE_H_
#ifdef WITH_ARRAY_TYPES
#ifdef HAVE_JAMA_JAMA_LU_H
#include "integrator.h"
#include "general.h"
#include "particle.h"
#include "MArray2D.h"
#include <string>
class Controller;
class Phase;
class IntegratorStaticLSE : public Integrator {

protected:
//
//	string  m_species2;
//
//	size_t m_colour2;
	  /*!
	   * The species of the two participating particles.
	   */
	  pair<string, string> m_species_pair;
	  /*!
	   * The name of the field
	   */
	  string m_scalar_name;
	  /*!
	   * The symbol (short name) of the field
	   */
	  string m_scalar_symbol;

	  string m_pairContribution_symbol;
	  /*!
	   * The symbol for the diagonal
	   */
	  string m_pairdiagonal_scalar_name;
	  /*!
	   * The tag offset of the scalar field
	   */
	  size_t m_pairdiagonal_scalar_offset;
	  /*!
	   * The tag offset of the scalar field
	   */
	  size_t m_particle_scalar_offset;

	  size_t m_pairContribution_symbol_offset;
	  /*!
	   *Boundary conditions
	   */
	  string m_boundary_condition_name;
	  /*
	   * The tag offset of the boundary condition
	   */
	  size_t m_boundary_condition_offset;

	  /*!
	   * The color pair belonging to the species combination.
	   */
      ColourPair *m_cp	;

      Pairdist* dis;
      ptrdiff_t m_pairslot;
//
//
      vector<PairList> m_freePairs;
      size_t nofixedParticles;

      /*! DofToMySlot tells us, which degree of freedom
       * of the system of equations corresponds to which slot
       */
	  size_t* DofToMySlot;

	  /*! mySlotToDof tells us, which particle slot
       * corresponds to which degree of freedom of the system of equations
       */
	  size_t* mySlotToDof;

 //matrix to be filled
	  MArray2D* m_matrixA;

	/*!
	 * Initialize the property list
	 */
 void init();
//	size_t getColourAndAdd(string species2);
public:
	/*!
	 * Constructor
	 * @param controller Pointer to the \a Controller object this \a Integrator belongs to
	 */
	IntegratorStaticLSE(Controller *controller);
	/*!
	 * Destructor
	 */
	virtual ~IntegratorStaticLSE();

	/*!
	 * Called right before the simulation will start
	 */

	/*!
	 * Add attributes
	 */
	virtual void setup();

	virtual void isAboutToStart();
	/*!
	 * Time integration step 1
	 */
	virtual void integrateStep1(){};
	/*!
	 * Time integration step 2
	 */
	virtual void integrateStep2();
	  /*!
	   * Return the tag offset of the scalar field
	   */
	  virtual size_t scalarOffset() const {
	    return m_particle_scalar_offset;
	  }
	  /*!
	   * Return the name of the scalar field
	   */
	  virtual string scalarName() const {
	    return m_scalar_name;
	  }

};
#endif /*WITH_JAMA_JAMA_LU*/
#endif /*HAVE_ARRAY_TYPES*/
#endif /*INTEGRATOR_STATIC_LSE_H_*/

