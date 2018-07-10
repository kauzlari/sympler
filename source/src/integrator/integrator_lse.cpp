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



#include "general.h"
#ifdef WITH_ARRAY_TYPES
#ifdef HAVE_JAMA_JAMA_LU_H
#include "integrator_lse.h"
#include "phase.h"
#include "controller.h"
#include "simulation.h"
#include "manager_cell.h"
#include "gen_f.h"
#include "threads.h"
#include <string>

#define CLASSNAME "IntegratorLSE"

#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()

const Integrator_Register<IntegratorLSE> integrator_lse("IntegratorLSE");

//---- Constructors/Destructor ----

IntegratorLSE::IntegratorLSE(Controller *controller) :
Integrator(controller) {
	init();
}

IntegratorLSE::~IntegratorLSE() {
}
//---- Methods ----
void IntegratorLSE::init() {
	m_properties.setClassName("IntegratorLSE");

	m_properties.setDescription("Integrates the linear system equations in time steps. "
			"The integrator assumes that dt is constant during the simulation.");

	DOUBLEPC(theta, m_theta, -INFINITY,
			"Theta parameter used in the trapezoidal rule for the integral calculation")
	;
	STRINGPC(symbol, m_symbol,
			"Input of the linear system equation")
	;
// 	STRINGPC(output_symbol, m_output_symbol,
// 			"Output of the linear system equation")
	;
	STRINGPC(terminalspecies, m_terminalspecies, "Terminal (e.g., wall) species")
	;
	BOOLPC(divideByDt, m_divideByDt, "Divide the right hand side by the timestep dt")
	;

	m_theta = 1;
	m_symbol = "UNDEF";
// 	m_output_symbol = "UNDEF";
	m_terminalspecies = "UNDEF";
	m_divideByDt = false; //was true
}

void IntegratorLSE::setup() {
	Integrator::setup();
	if (m_terminalspecies == "UNDEF") {
		throw gError("IntegratorLSE::setup",
				"No setup for terminalspecies UNDEF, please define terminalspecies.");
	}

	m_colour_lse = m_colour/*getColourAndAdd(m_terminalspecies)*/;
        m_colour = getColourAndAdd(m_terminalspecies);

	if (Particle::s_tag_format[m_colour_lse].attrExists("matrixA"))
		throw gError("pc_lse::init", "Attribute matrixA already exists");
	if (Particle::s_tag_format[m_colour_lse].attrExists("matrixE"))
		throw gError("pc_lse::init", "Attribute matrixE already exists");
	if (Particle::s_tag_format[m_colour_lse].attrExists("stateVector"))
		throw gError("pc_lse::init", "Attribute stateVector already exists");
	if (Particle::s_tag_format[m_colour_lse].attrExists("velocityVector"))
		throw gError("pc_lse::init", "Attribute velocityVector already exists");
	if (Particle::s_tag_format[m_colour].attrExists("matrixB"))
		throw gError("pc_lse::init", "Attribute matrixB already exist");
	if (Particle::s_tag_format[m_colour].attrExists("matrixC"))
		throw gError("pc_lse::init", "Attribute matrixC already exist");
	if (Particle::s_tag_format[m_colour].attrExists(m_symbol))
		throw gError("pc_lse::init", "Attribute "+ m_symbol +" already exists");
// 	if (Particle::s_tag_format[m_colour].attrExists(m_output_symbol))
// 		throw gError("pc_lse::init", "Attribute "+ m_output_symbol +" already exists");

	Particle::s_tag_format[m_colour_lse].addAttribute("matrixA",
			DataFormat::MArray2D, /*persist.first*/true, "matrixA");
	Particle::s_tag_format[m_colour_lse].addAttribute("matrixE",
			DataFormat::MArray2D, /*persist.first*/true, "matrixE");
	Particle::s_tag_format[m_colour_lse].addAttribute("stateVector",
			DataFormat::MArray2D, /*persist.first*/true, "stateVector");
	Particle::s_tag_format[m_colour_lse].addAttribute("velocityVector",
			DataFormat::MArray2D, /*persist.first*/true, "velocityVector");
	Particle::s_tag_format[m_colour].addAttribute("matrixB",
			DataFormat::MArray2D, /*persist.first*/true, "matrixB");
	Particle::s_tag_format[m_colour].addAttribute("matrixC",
			DataFormat::MArray2D, /*persist.first*/true, "matrixC");
	Particle::s_tag_format[m_colour].addAttribute(m_symbol,
			DataFormat::DOUBLE, /*persist.first*/true, m_symbol);
	for(size_t i = 0; i < FORCE_HIST_SIZE; ++i) {
	    Particle::s_tag_format[m_colour].addAttribute(
	    		STR_FORCE + STR_DELIMITER + m_symbol + STR_DELIMITER + ObjToString(i),
	    		DataFormat::DOUBLE, true);
	}
}

void IntegratorLSE::isAboutToStart() {
	MSG_DEBUG(string(CLASSNAME) << "::" << __func__,"Calculate inverse");

	// Only works when dt is constant during simulation!
	m_dt = M_CONTROLLER->dt();

	MSG_DEBUG(string(CLASSNAME) << "::" << __func__,
			"Species: " << m_species <<
			", terminalspecies: " << m_terminalspecies <<
			", Symbol: " << m_symbol/* <<
			", Output symbol: " << m_output_symbol*/ <<".");

	ParticleList* pl=&(M_PHASE->particles(m_colour_lse));
	offsetA=Particle::s_tag_format[m_colour_lse].attrByName("matrixA").offset;
	offsetE=Particle::s_tag_format[m_colour_lse].attrByName("matrixE").offset;
	offsetx=Particle::s_tag_format[m_colour_lse].attrByName("stateVector").offset;
	offsetv=Particle::s_tag_format[m_colour_lse].attrByName("velocityVector").offset;
	offsetB=Particle::s_tag_format[m_colour].attrByName("matrixB").offset;
	offsetC=Particle::s_tag_format[m_colour].attrByName("matrixC").offset;
	for(size_t i = 0; i < FORCE_HIST_SIZE; ++i) {

	    m_offset_force[i] = Particle::s_tag_format[m_colour].attrByName(
			"force_" + m_symbol + "_" + ObjToString(i)).offset;
	    m_fAttr_index[i] = Particle::s_tag_format[m_colour].attrByName(
			"force_" + m_symbol + "_" + ObjToString(i)).index;
        }
 /*Particle::s_tag_format[m_colour].attrByName(
			"force_" + m_input_symbol + "_" + ObjToString(i)).offset;*/
	m_offset_symbol=Particle::s_tag_format[m_colour].attrByName(m_symbol).offset;
	m_matrixA=*(*pl)[0].tag.array2dDoubleByOffset(offsetA);
	m_matrixE=*(*pl)[0].tag.array2dDoubleByOffset(offsetE);

	MSG_DEBUG(string(CLASSNAME) << "::" << __func__,
			"Dimensions: "
			"A(" << (m_matrixA.dim1()) << "," << (m_matrixA.dim2()) << "), "<<
			"E(" << (m_matrixE.dim1()) << "," << (m_matrixE.dim2()) << "); dt = "<<
			m_dt << " , theta = "<< m_theta <<".");

	if (m_theta != 0) {
		m_EA_inverted = invert((m_matrixE.scalarmult(1/(m_theta*m_dt)))
				-m_matrixA);
		//cout<<"the inverted matrix is"<<endl<<*m_EA_inverted.p_array2d;
	}
	else m_EA_inverted = invert(m_matrixE);

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (M_PHASE, m_colour, this,
       for (int j = 0; j < FORCE_HIST_SIZE; ++j)
         i->tag.doubleByOffset(((IntegratorLSE*) data)->m_offset_force[j]) = 0;
      );

}

void IntegratorLSE::integrateStep2() {
	size_t n = m_matrixA.dim1();
	MArray2D f(n, 1, 0);
	ParticleList* pl=&(M_PHASE->particles(m_colour_lse));
	MArray2D* x = ((*pl)[0].tag.array2dDoubleByOffset(offsetx)).value();
	MArray2D* v = ((*pl)[0].tag.array2dDoubleByOffset(offsetv)).value();
	size_t force_index = M_CONTROLLER->forceIndex();

	// Collect input
	FOR_EACH_FREE_PARTICLE_C__PARALLEL(M_PHASE,m_colour,0,
		for (int t = 0; t < global::n_threads; ++t) {
			f = f + (i->tag.array2dDoubleByOffset(offsetB)->scalarmult(
							i->tag.doubleByOffset(m_offset_force[force_index])));

// 			i->tag.doubleByOffset(m_offset_force[force_index]) = 0;


		}
	)
	;
	if (m_divideByDt)
		f = f.scalarmult(1/m_dt);
	//cout << (size_t)sbrk(0)<<endl;;

	// Integrate
	MArray2D xnew(n, 1);
	if(m_theta!=0) {
		xnew = matmult(m_EA_inverted,
				matmult(m_matrixE,
				((*x).scalarmult(1/m_theta/m_dt))+((*v).scalarmult((1/m_theta)-1)))-f);
		*v = (xnew-*x).scalarmult(1/m_theta/m_dt)-(*v).scalarmult(1/m_theta-1);
	}
	else {
		xnew = matmult(m_EA_inverted, matmult(m_matrixA,*x) + f).scalarmult(m_dt)+(*x);
		*v = (xnew-*x).scalarmult(1/m_dt);
	}

	// Assign output
	*x = xnew;
	FOR_EACH_FREE_PARTICLE_C__PARALLEL(M_PHASE,m_colour,0,
			i->tag.doubleByOffset(m_offset_symbol) =
			matmult(*(i->tag.array2dDoubleByOffset(offsetC)),*x)(0,0);
	)
	;
}


void IntegratorLSE::unprotect(size_t index)
{
  Phase *phase = M_PHASE;
  FOR_EACH_FREE_PARTICLE_C__PARALLEL
      (phase, m_colour, this,
       i->tag.unprotect(m_fAttr_index[index]);
       if(!index) i->tag.protect(m_fAttr_index[FORCE_HIST_SIZE-1]);
       else i->tag.protect(m_fAttr_index[index-1]);
      );

}


#ifdef _OPENMP
string IntegratorLSE::dofIntegr() {
  return m_symbol;
}


void IntegratorLSE::mergeCopies(Particle* p, int thread_no, int force_index) {
  if (m_merge == true) {
    p->tag.doubleByOffset(m_offset_force[force_index]) += (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos];
    (*p->tag.vectorDoubleByOffset(m_vec_offset[thread_no]))[m_vec_pos] = 0;
  }
}
#endif

#endif /*WITH_JAMA_JAMA_LU*/
#endif /*HAVE_ARRAY_TYPES*/
