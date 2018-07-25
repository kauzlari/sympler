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
#include "phase.h"
#include "controller.h"
#include "simulation.h"
#include "manager_cell.h"
#include "gen_f.h"
#include "threads.h"
#include <string>
#include"integrator_static_lse.h"

#ifdef WITH_INTEGRATORSTATICLSE
#ifdef WITH_ARRAY_TYPES
#ifdef HAVE_JAMA_JAMA_LU_H
#include "superlu/slu_ddefs.h"


#define CLASSNAME "IntegratorStaticLSE"

#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()

const Integrator_Register<IntegratorStaticLSE> integrator_static_lse(
		"IntegratorStaticLSE");

//---- Constructors/Destructor ----

IntegratorStaticLSE::IntegratorStaticLSE(Controller *controller) :
	IntegratorScalar(controller) {
	init();
}
IntegratorStaticLSE::~IntegratorStaticLSE() {
}
//---- Methods ----
void IntegratorStaticLSE::init() {
	m_properties.setClassName("IntegratorStaticLSE");

	m_properties.setDescription(
			"Integrates the linear equation system built from the "
				"entries in the colour pairs");

	STRINGPC
	(pairContribution, m_pairContribution_symbol,
			"Symbol assigned to the pair interactions, usable in algebraic expressions, for the off-diagonal elements")
	;
	STRINGPC
	(particleContribution, m_pairdiagonal_scalar_name,
			"Symbol assigned to the pair interactions, usable in algebraic expressions, for the diagonal elements")
	;

	STRINGPC(boundaryCondition,m_boundary_condition_name,"Boundary condition flag of the equations")
	;

}

void IntegratorStaticLSE::setup() {

  IntegratorScalar::setup();

  DataFormat::attribute_t tmpAttr;

  m_boundary_condition_offset
    = Particle::s_tag_format[m_colour].addAttribute
    (m_boundary_condition_name, DataFormat::DOUBLE, true, m_boundary_condition_name).offset;
  
  m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species), M_MANAGER->getColour(m_species));
  
  m_cp->setNeedPairs(true);

}

void IntegratorStaticLSE::isAboutToStart() {

  IntegratorScalar::isAboutToStart();
    
//	MSG_DEBUG(string(CLASSNAME) << "::" << __func__,
//			"Species: " << m_species <<
//			", Scalar: " << m_scalar_name <<
//			", Symbol: " << m_scalar_symbol <<
//			",BoundaryCondition: "<< m_boundary_condition_name <<".");

  m_pairdiagonal_scalar_offset = Particle::s_tag_format[m_colour].attrByName(
			m_pairdiagonal_scalar_name).offset;

}

void IntegratorStaticLSE::integrateStep2() {

	Phase *phase = M_PHASE;
	//	size_t realSize;
	//	realSize = 0;
//	m_matrixA=new MArray2D(3,5,0);
//		cout<< *m_matrixA;
//		MSG_DEBUG("IntegratorStaticLSE:isAboutToStart","Matrix dimensions are: "
					//<< (m_matrixA->dim1()) << "," << (m_matrixA->dim2()));
	//MSG_DEBUG("IntegratorStaticLSE::integrateStep2","n_threads=" <<global::n_threads);
	//Number of free particles
	size_t nofixedParticles = 0;
		size_t j = 0;


		FOR_EACH_FREE_PARTICLE_C
		  (phase, m_colour,
		     nofixedParticles++;
		   );
		//FIXME!: where is delete? Why not new once and then resizing?
		DofToMySlot = new size_t[nofixedParticles];
		mySlotToDof = new size_t[phase->returnNofPartC(m_colour)];
		FOR_EACH_FREE_PARTICLE_C
		  (phase, m_colour,
		   if ((__iSLFE->tag.doubleByOffset(m_boundary_condition_offset)<=0))
		     {
		       DofToMySlot[j] = __iSLFE->mySlot;
		       mySlotToDof[__iSLFE->mySlot] = j;
		       
		       ++j;
		     }
		   );

	FOR_EACH_PAIR(m_cp,

			DataFormat::attribute_t Attr =
			m_cp->tagFormat().attrByName(m_pairContribution_symbol);

	m_pairContribution_symbol_offset = Attr.offset; // TODO: For each colour pair separately
   //MSG_DEBUG("IntegratorStaticLSE::integrateStep2","PairContribution " <<(Attr.name));
  // MSG_DEBUG("IntegratorStaticLSE::integrateStep2","PairContribution " <<pair->tag.doubleByOffset(m_pairContribution_symbol_offset));
 // MSG_DEBUG("IntegratorStaticLSE::integrateStep2","FirstPartinPair " <<pair->firstPart()->mySlot);
  //MSG_DEBUG("IntegratorStaticLSE::integrateStep2","SecondPartinPair " <<pair->secondPart()->mySlot);
	);
	int numberOfParticles=M_PHASE->returnNofPartC(m_colour);
/*Initialization of the needed matrices*/

	if (nofixedParticles==numberOfParticles)
		m_matrixA= new MArray2D(nofixedParticles-1,nofixedParticles-1,0);

//	if(nofixedParticles != numberOfParticles)
//		m_matrixA= new MArray2D(nofixedParticles,nofixedParticles);
cout<<(*m_matrixA);

ParticleList* particles = &(M_PHASE->particles(m_colour));

//for all pairs, except the pairs in which the last particle is,fill the values
FOR_EACH_PAIR(m_cp,
	Particle* firstP = pair->firstPart();
    Particle* secondP = pair->secondPart();


    bool firstHasDirichletBC = firstP->tag.doubleByOffset(m_boundary_condition_offset)>0 || firstP -> isFrozen;
    bool secondHasDirichletBC = secondP->tag.doubleByOffset(m_boundary_condition_offset)>0 || secondP -> isFrozen;
    if ((!secondHasDirichletBC && !firstHasDirichletBC))
    {

if(mySlotToDof[pair->secondPart()->mySlot]<nofixedParticles-1 && mySlotToDof[pair->secondPart()->mySlot]<nofixedParticles-1)
{
    	(*m_matrixA)(pair->firstPart()->mySlot,pair->secondPart()->mySlot)= pair->tag.doubleByOffset(m_pairContribution_symbol_offset);
    	(*m_matrixA)(pair->secondPart()->mySlot,pair->firstPart()->mySlot)= pair->tag.doubleByOffset(m_pairContribution_symbol_offset);

    	//MSG_DEBUG("IntegratorStaticLSE::integrateStep2","index1 " << pair->firstPart()->mySlot);
    	//MSG_DEBUG("IntegratorStaticLSE::integrateStep2","index2 " << pair->secondPart()->mySlot);
    	//MSG_DEBUG("IntegratorStaticLSE:integrateStep2","eintrag"<<(*m_matrixA)(mySlotToDof[pair->firstPart()->mySlot],mySlotToDof[pair->secondPart()->mySlot]));
}
    }


);
//for all free particles but the last one,fill the diagonals
 FOR_EACH_FREE_PARTICLE_C
   (M_PHASE,m_colour,
    if(__iSLFE->tag.doubleByOffset(m_boundary_condition_offset)<=0 && __iSLFE->mySlot<nofixedParticles-1)
      {
	(*m_matrixA)(__iSLFE->mySlot,__iSLFE->mySlot)= __iSLFE->tag.doubleByOffset(m_pairdiagonal_scalar_offset);
      }
    );

//
//	/*Initialisation of the matrices needed*/
	SuperMatrix A,L,U,B;
	int *perm_c; /* column permutation vector */
	int *perm_r; /* row permutations from partial pivoting */
	int info;
	superlu_options_t options;
	SuperLUStat_t stat;

	double *rhs;



	MSG_DEBUG("IntegratorStaticLSE", "num of particles is" << numberOfParticles);

	size_t numberOfNonfixedPairs = 0;
	FOR_EACH_PAIR(m_cp,
			//cout << pair->firstPart()->tag.doubleByOffset(m_boundary_condition_offset) << "..." << pair->secondPart()->tag.doubleByOffset(m_boundary_condition_offset) << endl;
			if((pair->firstPart()->tag.doubleByOffset(m_boundary_condition_offset)<=0)&&(pair->secondPart()->tag.doubleByOffset(m_boundary_condition_offset)<=0))
		numberOfNonfixedPairs++;

	);
	// Pairs without boundary condition
	//size_t numberOfNonfixedPairs= (m_cp->freePairs().size()) +
		//(m_cp->frozenPairs().size());// nofixedParticles*(nofixedParticles-1)/2;

	// To store the entries of the sparse matrix
	size_t* nonZeroCoordinates;
	double* nonZeroValue;
	nonZeroCoordinates = (size_t*)malloc(sizeof(size_t)*2*(numberOfNonfixedPairs+nofixedParticles));
	nonZeroValue = (double*)malloc(sizeof(double)*(numberOfNonfixedPairs+nofixedParticles));
	//size_t nonZeroCoordinates[numberOfNonfixedPairs+nofixedParticles][2];
	//double nonZeroValue[numberOfNonfixedPairs+nofixedParticles];
	cout << "nofixedParticles " <<  nofixedParticles << endl;
	cout << "numberOfNonfixedPairs" << numberOfNonfixedPairs << endl;
	int currentNonzeroNumber=0;

	// Maximum DoF index of first and second parts of pairs
	size_t max_firstPart=0;
	size_t max_secondPart=0;

	// Doing the right hand side
	cout << "Alloc "<< (8*numberOfParticles) << " bytes" << endl;
	cout << "Alloc "<< (sizeof(double)*numberOfParticles) << " bytes" << endl;
	cout << "Alloc "<< (sizeof(double)*numberOfParticles) << " bytes" << endl;
	cout << "Alloc "<< (sizeof(double)*numberOfParticles) << " bytes" << endl;
	rhs = (double*)malloc(sizeof(double)*numberOfParticles);
	if ( !rhs ) ABORT("Malloc fails for rhs[].");
	for (int b = 0; b <numberOfParticles; ++b) rhs[b] = 0;
	cout<<"# no fixed"<<nofixedParticles<<endl;

	// Insert first diagonal entry into sparse notation
	size_t currentDof = 0;
	nonZeroCoordinates[0+1]=nonZeroCoordinates[0+0]=currentDof;
	nonZeroValue[0]=(*particles)[DofToMySlot[currentDof]].tag.doubleByOffset(m_pairdiagonal_scalar_offset);
	currentNonzeroNumber = 1;

	FOR_EACH_PAIR(m_cp,

			// Dirichlet boundary conditions:
			// If one of the particles has a BC, we modify the system of equations
			// Otherwise, we create a new off-diagonal entry in the sparse notation.
			Particle* firstP = pair->firstPart();
			Particle* secondP = pair->secondPart();
			bool firstHasDirichletBC = firstP->tag.doubleByOffset(m_boundary_condition_offset)>0 || firstP -> isFrozen;
			bool secondHasDirichletBC = secondP->tag.doubleByOffset(m_boundary_condition_offset)>0 || secondP -> isFrozen;

			//if the first is defined value,second is free substract from rhs
			if(firstHasDirichletBC && !secondHasDirichletBC) {
				rhs[mySlotToDof[secondP->mySlot]]-=(pair->tag.doubleByOffset(m_pairContribution_symbol_offset))*(firstP->tag.doubleByOffset(m_scalar_offset));
				//cout << "(" << (pair->secondPart()->mySlot) << "->" << (pair->tag.doubleByOffset(m_pair_scalar_symbol_offset)) << "*" << (pair->firstPart()->tag.doubleByOffset(m_scalar_offset)) << endl;
			}

			//if the second is defined value,substract from rhs
			if(secondHasDirichletBC && !firstHasDirichletBC) {
				rhs[mySlotToDof[firstP->mySlot]]-=pair->tag.doubleByOffset(m_pairContribution_symbol_offset)*secondP->tag.doubleByOffset(m_scalar_offset);
				//cout << "(" << (pair->firstPart()->mySlot) << "->" << (pair->tag.doubleByOffset(m_pair_scalar_symbol_offset)) << "*" << (pair->secondPart()->tag.doubleByOffset(m_scalar_offset)) << endl;
			}

			//if both are free in the matrix, means myslot_array_c only free particles
			if ((!secondHasDirichletBC && !firstHasDirichletBC))
			{

				nonZeroCoordinates[2*currentNonzeroNumber+0]=mySlotToDof[pair->firstPart()->mySlot];
				nonZeroCoordinates[2*currentNonzeroNumber+1]=mySlotToDof[pair->secondPart()->mySlot];
				//of a change is there then fill the value for the diagonals which are missed otherwise
				// When switching from (a,?) to (a+1,?), there is a good chance
				// to insert a diagonal matrix entry
				if(nonZeroCoordinates[2*currentNonzeroNumber+0]!=nonZeroCoordinates[2*(currentNonzeroNumber-1)+0])
				{
					nonZeroCoordinates[2*currentNonzeroNumber+1+0] = nonZeroCoordinates[2*currentNonzeroNumber+0];
					nonZeroCoordinates[2*currentNonzeroNumber+1+1] = nonZeroCoordinates[2*currentNonzeroNumber+1];
					nonZeroValue[currentNonzeroNumber+1] = nonZeroValue[currentNonzeroNumber];

					++currentDof;
					nonZeroCoordinates[2*currentNonzeroNumber+1]=nonZeroCoordinates[2*currentNonzeroNumber+0]=currentDof;
					nonZeroValue[currentNonzeroNumber]=(*particles)[DofToMySlot[currentDof]].tag.doubleByOffset(m_pairdiagonal_scalar_offset);
					++currentNonzeroNumber;
				}
				nonZeroValue[currentNonzeroNumber]=pair->tag.doubleByOffset(m_pairContribution_symbol_offset);
//				cout << "(myslot_array_c[j],myslot_array_v[j]) = ((" << nonZeroCoordinates[currentNonzeroNumber][0] << "," <<
//				nonZeroCoordinates[currentNonzeroNumber][1] << ")," << nonZeroValue[currentNonzeroNumber] << ")" << endl;

				if (nonZeroCoordinates[2*currentNonzeroNumber+0]>max_firstPart)
				max_firstPart= nonZeroCoordinates[2*currentNonzeroNumber+0];

				if (nonZeroCoordinates[2*currentNonzeroNumber+1]>max_secondPart)
				max_secondPart=nonZeroCoordinates[2*currentNonzeroNumber+1];

				currentNonzeroNumber++;
			}
	);
	++currentDof;
	nonZeroCoordinates[2*currentNonzeroNumber+1]=nonZeroCoordinates[2*currentNonzeroNumber+0]=currentDof;
	nonZeroValue[currentNonzeroNumber]=(*particles)[DofToMySlot[currentDof]].tag.doubleByOffset(m_pairdiagonal_scalar_offset);

	int maxDofIndex=0;
	if (max_firstPart>max_secondPart)
	maxDofIndex=max_firstPart+1;
	else maxDofIndex=max_secondPart+1;

	//Convert to strange superlu format
	size_t numReallyFreePairs=currentNonzeroNumber+1; // +1 because of last diagonal element
	vector<int> nzr (maxDofIndex,0);

	for(int i=0;i<numReallyFreePairs;i++)
	{
		if (nonZeroValue[i]!= 0)
		{
			nzr.at(nonZeroCoordinates[2*i+0])++;

			if(nonZeroCoordinates[2*i+0]!=nonZeroCoordinates[2*i+1])
			{
				nzr[nonZeroCoordinates[2*i+1]]++;
			}
			else {}
		}
		else {}
	}

	//summing up to get the total number of non zeros in the array (nz)
	int nz=0;
	for (vector<int>::const_iterator itr = nzr.begin();
			itr != nzr.end();
			itr++)
	{
		nz += *itr;
	}
//	cout<<"Total number of non zeros in the array"<<nz<<endl;
	//	//xa - array for the indices indicating the beginning of each column in the coefficient and row index array
	//int xa[num1+1];
	int xa[maxDofIndex+1];
	xa[0]= 0;

	for(int j=1;j<maxDofIndex+1;j++)
	{
		xa[j]=xa[j-1]+nzr.at(j-1);
	}

	//asub[] - stores the row indices of the non-zero coefficients
	// a[] - stores the non zero coefficients of the array column numbered
	// ptr- intermediate vector
	int asub[nz];
	double a[nz];

	vector<int> ptr(maxDofIndex,1);

	for(size_t j=0;j<numReallyFreePairs;j++)
	{
		if (nonZeroValue[j]!= 0)
		{
			int s1=xa[(nonZeroCoordinates[2*j+0])]-1+ptr.at(nonZeroCoordinates[2*j+0]);
			a[s1]=nonZeroValue[j];
			//			cout<<"index0"<<myslot_array[j][0]<<" ptr"<<ptr.at(myslot_array[j][0])<<endl;
			//			cout<<"j is"<<j<<" s1 is"<<s1 <<" a value is"<<a[s1]<<endl;
			asub[s1]=nonZeroCoordinates[2*j+1];
			ptr.at(nonZeroCoordinates[2*j+0])++;

			if(nonZeroCoordinates[2*j+0]!=nonZeroCoordinates[2*j+1])
			{
				int s2=xa[(nonZeroCoordinates[2*j+1])]-1+ptr.at(nonZeroCoordinates[2*j+1]);
				a[s2]=nonZeroValue[j];

				asub[s2]=nonZeroCoordinates[2*j+0];
				ptr.at(nonZeroCoordinates[2*j+1])++;

			}
			else {}
		}
		else {}
	}

	//		 Create matrix A in the format expected by SuperLU.
	//dCreate_CompCol_Matrix(&A, nofixedParticles, nofixedParticles, nz, &a[0], &asub[0], &xa[0], SLU_NC, SLU_D, SLU_GE);
	//		 Create right-hand side matrix B.

	MSG_DEBUG("IntegratorStaticLSE", "num of particles is" << numberOfParticles);
//	for (size_t i=0;i<nofixedParticles;i++)
//	{
//		cout<<"rechte seite:"<<rhs[i]<<endl;
//	}
	//dCreate_Dense_Matrix(&B,nofixedParticles, 1, rhs, nofixedParticles, SLU_DN, SLU_D, SLU_GE);
	if ( !(perm_r = intMalloc(nofixedParticles))) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(nofixedParticles))) ABORT("Malloc fails for perm_c[].");
	//	 Set the default input options.


	set_default_options(&options);
	options.ColPerm = NATURAL;
	//		//	 Initialize the statistics variables.
	StatInit(&stat);
	MSG_DEBUG("IntegratorStaticLSE::integrateStep2","before calculation");
	//dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
	MSG_DEBUG("IntegratorStaticLSE::integrateStep2","after calculation");
	//dPrint_CompCol_Matrix("A", &A);

	//dPrint_CompCol_Matrix("U", &U);
	//dPrint_SuperNode_Matrix("L", &L);
	//print_int_vec("\nperm_r", m, perm_r);
	//dPrint_Dense_Matrix("The solution matrix B is",&B);
	//cout<<"number of rows in B  "<<B.nrow<<endl;
	DNformat *Astore;
	double *solution_vec;
	//Astore = (DNformat *) B.Store;
	solution_vec = (double *) Astore->nzval;

	int s=0;

	// Store solution in particle DoF
	FOR_EACH_FREE_PARTICLE_C
	  (M_PHASE,m_colour,
	   if(__iSLFE->tag.doubleByOffset(m_boundary_condition_offset)<=0)
	     {
	       __iSLFE->tag.doubleByOffset(m_scalar_offset)=solution_vec[s];
	       s++;
	     }
	   );

	// FIXME!: question to the main authors cenovai and lieneman: 
	// what do the next lines mean? TODO or NOT TODO?
	//	/* TODO: De-allocate storage */To be included when completed
	//	SUPERLU_FREE (rhs);
	//	SUPERLU_FREE (perm_r);
	//	SUPERLU_FREE (perm_c);
	//	Destroy_CompCol_Matrix(&A);
	//	Destroy_SuperMatrix_Store(&B);
	//	Destroy_SuperNode_Matrix(&L);
	//	Destroy_CompCol_Matrix(&U);
	//	StatFree(&stat);
}


#endif /*WITH_JAMA_JAMA_LU*/
#endif /*HAVE_ARRAY_TYPES*/
#endif /*IntegratorStaticLSE*/

