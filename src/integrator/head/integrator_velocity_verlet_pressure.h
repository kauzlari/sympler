/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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



#ifndef INTEGRATOR_VELOCITY_VERLET_PRESSURE_H_
#define INTEGRATOR_VELOCITY_VERLET_PRESSURE_H_

#ifdef WITH_ARRAY_TYPES

#include "integrator_velocity_verlet.h"
#include "particle.h"
#include "MArray2D.h"
#include <string>
#include "particle.h"
#include "MArray2D.h"
#include <string>
#include "general.h"

using namespace std;

class GenF;
class Phase;
class Controller;
class WallTriangle;
class Cell;

//----IntegratorVelocityVerlet ----

/*!
 * Modified Velocity-Verlet integrator for the positions and velocities
 * See: R. D. Groot and P. B. Warren, J. Chem. Phys. 107, 4423-4435 (1997)
 */

class IntegratorVelocityVerletPressure: public IntegratorVelocityVerlet {
protected:

	/*Symbol for nabla of the weighting function*/
	string m_nablaWF_symbol;
	/*!
	 * The symbol for the scalar diagonal in the LAPLACE coefficient matrix
	 */
	//string m_LaplaceParticle_name;
	/*!
	 * The tag offset of the scalar field
	 */
	size_t m_nablaWF_offset;
	/*!
	 * The species of the two participating particles.
	 */
	pair<string, string> m_species_pair;
	/*!
	 * The name of the scalar field
	 */
	string m_vector_name;

	/*!
	 * The symbol (short name) of the scalar field
	 */
	string m_vector_symbol;

	/*!
	 * The tag offset of the scalar field
	 */
	size_t m_vector_offset;
	/*!
	 * The name of the scalar field
	 */
	string m_velcorr_name;

	/*!
	 * The symbol (short name) of the scalar field
	 */
	string m_velcorr_symbol;

	/*!
	 * The tag offset of the scalar field
	 */
	size_t m_velcorr_offset;

	double rho;
	string m_rho_name;
	size_t m_rho_offset;

	size_t m_tensor_offset;
	/*!
	 * Symbol of the given tensor
	 */
	string m_tensor_symbol;
	/*!
	 * Initialize the property list
	 */
	void init();
	/*!
	 * The color pair belonging to the species combination.
	 */
	ColourPair *m_cp;

	Pairdist* dis;
	ptrdiff_t m_pairslot;
	//matrix to be filled
	MArray2D m_matrixD_ijx;
	MArray2D m_matrixD_ijy;
	MArray2D m_matrixD_ijz;

	MArray2D m_matrixR_ijx;
	MArray2D m_matrixR_ijy;
	MArray2D m_matrixR_ijz;

	MArray2D m_matrixL_ikx;
	MArray2D m_matrixL_iky;
	MArray2D m_matrixL_ikz;

	MArray2D m_matrixL;

	MArray2D grad;
	MArray2D* result_c;
	int numberOfParticles;
	MArray2D div_v;
	MArray2D result_p;
	MArray2D grad_v;

public:
	/*!
	 * Constructor
	 * @param controller Pointer to the \a Controller object this \a Integrator belongs to
	 */
	IntegratorVelocityVerletPressure(Controller *controller);

	/*!
	 * Destructor
	 */
	virtual ~IntegratorVelocityVerletPressure();
	/*!
	 * Add attributes
	 */
	virtual void setup();
	/*!
	 * Initialize temporary fields and clear all forces
	 */
	virtual void isAboutToStart();

	/*!
	 * Estimation of the velocity
	 */
	virtual void integrateStep1();

	/*!
	 * Correction of the velocity
	 */
	virtual void integrateStep2();

	/*!
	 * Integration of the position
	 */
	virtual void integratePosition(Particle* p, Cell* cell);

	/*!
	 * Prediction of the velocity
	 */
	virtual void integrateVelocity(Particle* p);

	/*!
	 * Solves the equation that checks for hits
	 */
	virtual void solveHitTimeEquation(WallTriangle* wallTriangle,
			const Particle* p, const point_t &force, vector<double>* results);

	/*!
	 * Checks which of the times (in the time vector) is the actual hit position. The function will
	 be used in WallTriangle
	 */
	virtual void hitPos
	  (const double& dt, const Particle* p, point_t &hit_pos, const point_t &force);

#ifdef _OPENMP

	// all here inherited from parent 
	
#endif

};

#endif /* WITH_ARRAY_TYPES */

#endif /* INTEGRATOR_VELOCITY_VERLET_PRESSURE_H_ */
