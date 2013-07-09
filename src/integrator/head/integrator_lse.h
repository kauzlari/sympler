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



#ifndef INTEGRATOR_LSE_H_
#define INTEGRATOR_LSE_H_
#include "general.h"
#ifdef WITH_ARRAY_TYPES
#ifdef HAVE_JAMA_JAMA_LU_H
#include "integrator.h"
#include "particle.h"
#include "MArray2D.h"
#include <string>
class Controller;
class Phase;
class IntegratorLSE : public Integrator {
protected:
	size_t m_colour_lse;

        /*!
         * The indices of the force attribute of the tag of the force on this scalar field.    */
        size_t m_fAttr_index[FORCE_HIST_SIZE];

	MArray2D m_matrixA;
	MArray2D m_matrixE;
	MArray2D m_EA_inverted;
	/*!
	 * Time step
	 */
	double m_dt, m_theta;
	string m_symbol, /*m_output_symbol*/ m_terminalspecies;
	size_t offsetA, offsetE, offsetx, offsetv, m_offset_symbol, offsetB,
			offsetC;
	size_t m_offset_force[FORCE_HIST_SIZE];

	bool m_divideByDt;
public:
	/*!
	 * Constructor
	 * @param controller Pointer to the \a Controller object this \a Integrator belongs to
	 */
	IntegratorLSE(Controller *controller);
	/*!
	 * Destructor
	 */
	virtual ~IntegratorLSE();
	/*!
	 * Initialize the property list
	 */
	void init();

        /*!
         * If used forces are saved in a tag, protect and unprotect them as needed
         */
        virtual void unprotect(size_t index);

#ifdef _OPENMP
       /*!
        * Returns the number of doubles this Integrator saves in a particle tag.
        */
        virtual int numCopyDoubles() {
          return 1;
        }

        virtual string dofIntegr();

        /*!
        * Merge the copies at the end of every timestep
        */
        virtual void mergeCopies(Particle* p, int thread_no, int force_index);

#endif

	/*!
	 * Called right before the simulation will start
	 */
	virtual void isAboutToStart();
	/*!
	 * Time integration step 1
	 */
	virtual void integrateStep1() {
	}
	/*!
	 * Time integration step 2
	 */
	virtual void integrateStep2();

	/*!
	 * Add attributes
	 */
	virtual void setup();
};
#endif /*WITH_JAMA_JAMA_LU*/
#endif /*HAVE_ARRAY_TYPES*/
#endif /*INTEGRATOR_LSE_H_*/
