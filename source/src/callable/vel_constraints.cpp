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



#include "vel_constraints.h"

#include "phase.h"
#include "random.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"
#ifdef HAVE_JAMA_JAMA_LU_H
#include "jama/jama_lu.h"
using namespace std;

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

const Callable_Register<VelConstraints> vel_constraints("VelConstraints");
/* ---VelConstraints --- */

VelConstraints::VelConstraints(Simulation* sim) :
	Callable(sim) {
	init();
}

VelConstraints::~VelConstraints() {
}

void VelConstraints::init() {
	m_properties.setClassName("VelConstraints");

	m_properties.setDescription("Velocity imposing constraints.");

	STRINGPC
	(species, m_species,
			"Species this VelConstraints should work on.")
		;

	STRINGPC
	(expression, m_exprString,
			"Algebraic vector expression for the position of the constrained particles.")
		;
	STRINGPC
	(vels_expression, m_vels_exprString,
			"Algebraic vector expression for the constrained velocity.")
		;
	m_species = "UNDEF";
	m_exprString = "idVec(1)";
	STRINGPC
	(symbol, m_constr_symbol,
			"Symbol assigned to the gradient P, usable in algebraic expressions")
		;

	STRINGPC
	(velcorr, m_constr_name,
			"Full name of the additional scalar field, usable as attribute in other modules")
		;
}

void VelConstraints::setup() {
	Callable::setup();

	if (m_species == "UNDEF")
		throw gError("VelConstraints::setup",
				"Attribute 'species' was not defined!");

	m_colour = M_MANAGER->getColour(m_species);

	m_expression.setExpression(m_exprString);
	m_expression.setReturnType(Variant::VECTOR);
	m_expression.setColour/*Pair*/(/*m_cp,*/m_colour);

	m_vels_expression.setExpression(m_vels_exprString);
	m_vels_expression.setReturnType(Variant::VECTOR);
	m_vels_expression.setColour/*Pair*/(/*m_cp,*/m_colour);
	//MSG_DEBUG("VelConstraints::setup", "expression= " << expression);
	m_constr_offset = Particle::s_tag_format[m_colour].addAttribute(
			m_constr_name, DataFormat::POINT, true, m_constr_symbol).offset;
}
void VelConstraints::call(size_t timestep) {
	constrain(M_PHASE);
}
void VelConstraints::constrain(Phase* phase) {
	velCMc[0] = 0;
	velCMc[1] = 0;
	velCMc[2] = 0;

	point_t temp;
	point_t vels_temp;
	int partnum = 0;
	FOR_EACH_FREE_PARTICLE_C
	  (phase,m_colour,
	   m_expression(&temp, __iSLFE);
	   if (temp[0] !=0 || temp[1] !=0 || temp[2] !=0)
	     {
	       velCMc += __iSLFE->v;
	       partnum++;
	     }
	   );

	MArray2D corr(partnum, 1, 0);
	point_t idvec;
	idvec[0] = 1;
	idvec[1] = 0;
	idvec[2] = 0;
	double value = 0;

	//MSG_DEBUG("VelConstraints::constrain", "constrpart " << partnum);
	//MSG_DEBUG("VelConstraints::constrain", "expression " << temp);
	for (int j = 0; j < SPACE_DIMS; j++) {
		velCMc[j] /= partnum;
	}
	//MSG_DEBUG("VelConstraints::constrain", "partnum" << partnum);
	// apply needed correction to desired velocity
	FOR_EACH_FREE_PARTICLE_C
	  (phase,	m_colour,
	   m_vels_expression(&vels_temp, __iSLFE);
	   m_expression(&temp, __iSLFE);
	   
	   if (temp[0] !=0 || temp[1] !=0 || temp[2] !=0)
	     __iSLFE->v += vels_temp-velCMc;
	   value+=(vels_temp[2]-velCMc[2]);
	   
	   __iSLFE->tag.pointByOffset (m_constr_offset)=vels_temp-velCMc;
	   );

	value = value / partnum;

}
#endif
