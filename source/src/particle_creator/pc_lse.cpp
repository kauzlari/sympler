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
#include "pc_lse.h"

#define CLASSNAME "ParticleCreatorLSE"

#define M_BOUNDARY ((Boundary*) m_parent)
#define M_PHASE ((Phase*) M_BOUNDARY->parent())
#define M_SIMULATION ((Simulation*) M_PHASE->parent())
#define M_MANAGER M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

/* Register this ParticleCreator with the factory. */
const ParticleCreator_Register<ParticleCreatorLSE>
		particle_creator_LSE("ParticleCreatorLSE");

//constructor
ParticleCreatorLSE::ParticleCreatorLSE(Boundary *boundary) :
	ParticleCreator(boundary) {
	init();
}
//destructor
ParticleCreatorLSE::~ParticleCreatorLSE() {

}

void ParticleCreatorLSE::setup() {
	ParticleCreator::setup();
	MSG_DEBUG(string(CLASSNAME)+"::"+__func__,
			"Species: "+m_species+", Terminalspecies: "+m_terminalspecies);
	//size_t terminalcolour=M_MANAGER->getColour(m_terminalspecies);
	//size_t m_colour=M_MANAGER->getColour(m_species);
}
//Methods

void ParticleCreatorLSE::createParticles() {
	//Controller* c = M_CONTROLLER;
	MSG_DEBUG(CLASSNAME+("::"+PRETTYFUNC_INFO),"Read matrices");
	M_CONTROLLER->registerForSetupAfterParticleCreation(this);
	if (m_matfile_A == "")
		throw gError("ParticleCreatorLSE::createParticles", "Please define the name of the matrix file A.");
	if (m_matfile_B== "")
		throw gError("ParticleCreatorLSE::createParticles", "Please define the name of the matrix file B.");
	if (m_matfile_C == "")
		throw gError("ParticleCreatorLSE::createParticles", "Please define the name of the matrix file C.");
	if (m_matfile_E == "")
		throw gError("ParticleCreatorLSE::createParticles", "Please define the name of the matrix file E.");
	m_matrixA=new MArray2D(m_matfile_A);
	m_matrixB=new MArray2D(m_matfile_B);
	m_matrixC=new MArray2D(m_matfile_C);
	m_matrixE=new MArray2D(m_matfile_E);
	*m_matrixA = (m_matrixA)->scalarmult(m_mult_matA);
	*m_matrixB = (m_matrixB)->scalarmult(m_mult_matB);
	*m_matrixC = (m_matrixC)->scalarmult(m_mult_matC);
	*m_matrixE = (m_matrixE)->scalarmult(m_mult_matE);
	MSG_DEBUG(CLASSNAME+("::"+PRETTYFUNC_INFO),"Matrix dimensions are: "
			"A(" << (m_matrixA->dim1()) << "," << (m_matrixA->dim2()) << "), "<<
			"E(" << (m_matrixE->dim1()) << "," << (m_matrixE->dim2()) << "), "<<
			"B(" << (m_matrixB->dim1()) << "," << (m_matrixB->dim2()) << "), "<<
			"C(" << (m_matrixC->dim1()) << "," << (m_matrixC->dim2()) << ")."
	);
	MSG_DEBUG(CLASSNAME+("::"+PRETTYFUNC_INFO),"Create superparticle");
	Particle p;
	p.setColour(m_colour);
	size_t offsetA=Particle::s_tag_format[m_colour].attrByName("matrixA").offset;
	size_t offsetE=Particle::s_tag_format[m_colour].attrByName("matrixE").offset;
	size_t offsetx=Particle::s_tag_format[m_colour].attrByName("stateVector").offset;
	size_t offsetv=Particle::s_tag_format[m_colour].attrByName("velocityVector").offset;
	p.tag.array2dDoubleByOffset(offsetA).set(m_matrixA);
	p.tag.array2dDoubleByOffset(offsetE).set(m_matrixE);
	p.tag.array2dDoubleByOffset(offsetx).set(new MArray2D(m_matrixA->dim1(),1,0));
	p.tag.array2dDoubleByOffset(offsetv).set(new MArray2D(m_matrixA->dim1(),1,0));
	M_PHASE->addParticle(p);
	MSG_DEBUG(CLASSNAME+("::"+PRETTYFUNC_INFO),"Done");

}
//Each terminal species gets a column form B matrix and a row form the C matrix.
void ParticleCreatorLSE::setupAfterParticleCreation() {
	size_t colour=M_MANAGER->getColour(m_terminalspecies);
	//ParticleList* p=&(M_PHASE->particles(c));
	//makro statt loop,the list durchgehen und for each particle attribute hinzufügen (zeile von C, spalte von B)
	//FOR_EACH loop over a vector,list,ect
	size_t offsetB, offsetC;
	offsetB=Particle::s_tag_format[colour].attrByName("matrixB").offset;
	offsetC=Particle::s_tag_format[colour].attrByName("matrixC").offset;
	int j=0;

	MSG_DEBUG(CLASSNAME+("::"+PRETTYFUNC_INFO),"Assign B and C matrices to particles");
	FOR_EACH_FREE_PARTICLE_C__PARALLEL(M_PHASE,colour,0,
			i->tag.array2dDoubleByOffset(offsetB) =
			array2d_double_sp(&(m_matrixB->column(j)));
			i->tag.array2dDoubleByOffset(offsetC) =
			array2d_double_sp(&(m_matrixC->row(j)));
			j++;
	)
	;
	//			FOR_EACH_FREE_PARTICLE_C__PARALLEL(M_PHASE,colour,0,
	//			cout<<"the added matrix element B is"<<*(i->tag.array2dDoubleByOffset(offsetB))
	//			;
	//			cout<<"the added matrix element C is"
	//					<<*(i->tag.array2dDoubleByOffset(offsetC));
	//			j++;
	//			)
	//			;
}

void ParticleCreatorLSE::init() {
	//ParticleCreator::init();
	m_properties.setClassName("ParticleCreatorLSE");
	m_properties.setDescription("Generates a super particle");
	STRINGPC(matrixA, m_matfile_A, "Coefficient matrix A")
	;
// 	default
	m_matfile_A="";
	STRINGPC(matrixE, m_matfile_E, "Free coefficient matrix E")
	;
	m_matfile_E="";
	STRINGPC(matrixB, m_matfile_B, "Coefficient matrix B")
	;
	m_matfile_B="";
	STRINGPC(matrixC, m_matfile_C, "Free coefficient matrix C")
	;
	m_matfile_C="";
	STRINGPC(terminalspecies, m_terminalspecies, "Terminal (e.g., wall) species")
	;

	m_mult_matA = 1;
	DOUBLEPC(multMatrixA, m_mult_matA, -INFINITY, "Multiplier of matrix A")
	;
	m_mult_matE = 1;
	DOUBLEPC(multMatrixE, m_mult_matE, -INFINITY, "Multiplier of matrix E")
	;
	m_mult_matB = 1;
	DOUBLEPC(multMatrixB, m_mult_matB, -INFINITY, "Multiplier of matrix B")
	;
	m_mult_matC = 1;
	DOUBLEPC(multMatrixC, m_mult_matC, -INFINITY, "Multiplier of matrix C")
	;
}

ostream &ParticleCreatorLSE::write(ostream &s, int shift) {
	return Node::write(s, shift);
}
#endif //WITH_ARRAY_TYPES
