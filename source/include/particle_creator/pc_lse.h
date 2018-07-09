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



#ifndef PC_LSE_H_
#define PC_LSE_H_
#include "general.h"
#ifdef WITH_ARRAY_TYPES
#include "phase.h"
#include "node.h"
#include "boundary.h"
#include "manager_cell.h"
#include "particle_creator.h"
#include "threads.h"
#include "simulation.h"

/* ---- ParticleCreatorLSE ---- */

/*!
 * Generates a super particle.
 */
class ParticleCreatorLSE : public ParticleCreator {
protected:

	/*!
	 * Initialisation of the \a PropertyList
	 */
	void init();

	/*!
	 * Setup for the \a ParticleCreatorLSE
	 */
	virtual void setup();
	//name of the files to be read
	string m_matfile_A;
	string m_matfile_B;
	string m_matfile_C;
	string m_matfile_E;

	// multipliers for matrices
	double m_mult_matA;
	double m_mult_matB;
	double m_mult_matC;
	double m_mult_matE;

	string m_terminalspecies;
	//name of the m_property
	MArray2D* m_matrixA;
	MArray2D* m_matrixB;
	MArray2D* m_matrixC;
	MArray2D* m_matrixE;
	MArray2D* m_matrix;

public:
	/*!
	 * Constructor
	 * @param boundary The \a Boundary, this \a ParticleCreator belongs to
	 */
	ParticleCreatorLSE(Boundary *boundary);
	/*!
	 * Destructor
	 */
	virtual ~ParticleCreatorLSE();

	/*!
	 * The routine for creating \a Particle s
	 */

	virtual void createParticles();
	virtual void setupAfterParticleCreation();

	/*!
	 * Write out the final state
	 */
	virtual ostream& write(ostream &s, int shift);
};

#endif /*WITH_ARRAY_TYPES*/
#endif /*PC_LSE_H_*/
