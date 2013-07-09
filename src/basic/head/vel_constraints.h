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



#ifndef __VEL_CONSTRAINTS_H
#define __VEL_CONSTRAINTS_H

#include "random.h"
#include "general.h"
#include "callable.h"
#include "colour_pair.h"
#include "function_particle.h"
#include "MArray2D.h"

using namespace std;

class Phase;
class Simulation;

/*!
 * A \a VelConstraints is a class for velocity constraints, i.e., modules that keep the mean velocity of
 * the particles to one defined value
 */
class VelConstraints: public Callable {
protected:

	/*!
	 * The colour, this \a VelConstraints should act on
	 */
	size_t m_colour;
	/*!
	 * The species, this \a VelConstraints should act on
	 */
	string m_species;
	/*!
	 * The center of mass velocity
	 */
	point_t velCMc;

	FunctionParticle m_expression;
	/*!
	 * The string for \a m_expression
	 */
	string m_exprString;

	FunctionParticle m_vels_expression;
	/*!
	 * The string for \a m_expression
	 */
	string m_vels_exprString;
	/*!
	 * The name of the scalar field
	 */
	string m_constr_name;

	/*!
	 * The symbol (short name) of the scalar field
	 */
	string m_constr_symbol;
	/*!
	 * The tag offset of the scalar field
	 */
	size_t m_constr_offset;
	/*!
	 * Initialize the property list
	 */
	void init();

public:
	/*!
	 * Constructor
	 * @param sim Pointer to the main simulation object
	 */
	VelConstraints(Simulation* sim);

	/*!
	 * Destructor
	 */
	virtual ~VelConstraints();

	/*!
	 * This function calls the velocity constraint, do I need it????
	 */
	virtual void call(size_t timestep);

	/*!
	 * Setup this \a VelConstraints
	 */
	virtual void setup();
	/*!
	 * Constrain the system
	 */
	virtual void constrain(Phase* p);
};

#endif /* __VEL_CONSTRAINTS_H */
