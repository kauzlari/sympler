/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2016, 
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

#include <stdlib.h>

#include <ostream>
#include <algorithm>

#include "simulation.h"
#include "manager_cell.h"
#include "threads.h"

#include "lin_solve_tensor_vector.h"

using namespace std;

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

const Callable_Register<LinSolveTensorVector> lin_solve_tensor_vector("LinSolveTensorVector");

LinSolveTensorVector::LinSolveTensorVector(Simulation* sim)
  : Callable(sim), m_cps(3,NULL), m_functions(3,NULL), m_1stParticleFactors(3,NULL), m_2ndParticleFactors(3,NULL), m_mat(), m_countdown(0), m_step_counter(0)
{
  init();
}

LinSolveTensorVector::~LinSolveTensorVector()
{
}

void LinSolveTensorVector::init()
{
  m_properties.setClassName("LinSolveTensorVector");
  
  m_properties.setDescription
    ("A LinSolveTensorVector solves the linear system of equations K.v=F.\n" 
     "* User defined input: tensor pair symbol for off-diagonal terms of matrix K, name for vector pair symbol holding solution v, vector pair symbol for right-hand side F. The diagonal entry in a row of K is the negative sum of the off-diagonal terms of this row. If species1 != species2, the module creates a linear system over all free particles of both species.\n"
     "* Optional output: matrix K in ASCII or BINARY format.\n"
     "* FIXME: 1) The module does not yet solve but just outputs the matrix\n"
     "* FIXME: 2) Matrix is currently stored non-sparse. So don't use with more than 3000 particles!");

  STRINGPC
    (species1, m_species1, 
     "First particle species."
     );
  m_species1 = "undefined";

  STRINGPC
    (species2, m_species2, 
     "Second particle species."
     );
  m_species2 = "undefined";

  STRINGPC
    (solSymbolName, m_solSymbolName, 
     "Name of the symbol storing the solution of the linear system. Symbol will be created by this module."
     );
  m_solSymbolName = "undefined";

  STRINGPC
      (pairFactor, m_expression,
       "The mathematical expression to be computed for the pairFactor."
      );
  m_expression = "undefined";
  
  STRINGPC
      (particleFactor_i, m_1stPExpression,
       "The mathematical expression of the additional particle factor for the first particle."
      );
  m_1stPExpression = "undefined";

  STRINGPC
      (particleFactor_j, m_2ndPExpression,
       "The mathematical expression of the additional particle factor for the second particle."
      );
  m_2ndPExpression = "undefined";

  STRINGPC
      (RHS, m_RHSsymbolName,
       "Name of the symbol holding the right hand side of the equation. Symbol must already exist."
      );
  m_RHSsymbolName = "undefined";

  INTPC
    (symmetry, m_symmetry, -2,
     "Factor of 1 or -1 that makes the matrix symmetric or antisymmetric"
     );
  m_symmetry = 0;

  INTPC
    (writeSystemEvery, m_writeSystemEvery, -1,
       "How often to write the matrix to a file. 0 means never, 1 in every time step."
      );
  m_writeSystemEvery = 0;

  DOUBLEPC
    (cutoff, m_cutoff, 0,
       "Cutoff for the pair interaction. Note that it should be consistent with previously computed quantities entering one of the runtime compiled expressions of this module. Must be >0."
      );
  m_writeSystemEvery = 0;

  BOOLPC
    (binary, m_binary,
     "If set to \"no\", \"false\", \"0\", ascii files are produced, if set to \"yes\", \"true\", \"1\", binary files. In binary mode, of course no header is written. Is ignored if m_writeSystemEvery = \"0\"."
     );
  m_binary = 0;

  STRINGPC
      (nameOutputFile, m_outputFileName,
       "Name of file for writing the matrix. The RHS is written into a file with postfix \"RHS\". Is ignored if writeSystemEvery = \"0\"."
      );
  m_outputFileName = "undefined";

} 

void LinSolveTensorVector::setup()
{

  if(m_species1 == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'species1' has value \"undefined\".");

  if(m_species2 == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'species2' has value \"undefined\".");

  if(m_solSymbolName == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'solSymbolName' has value \"undefined\".");

  if(m_expression == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'pairFactor' has value \"undefined\".");

  if(m_1stPExpression == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'particleFactor_i' has value \"undefined\".");

  if(m_2ndPExpression == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'particleFactor_j' has value \"undefined\".");

  if(m_RHSsymbolName == "undefined")
    throw gError("LinSolveTensorVector::setup", "Attribute 'RHS' has value \"undefined\".");

  if(m_symmetry != 1 && m_symmetry != -1)
    throw gError("LinSolveTensorVector::setup", "Attribute 'symmetry' has value \"" + ObjToString(m_symmetry) + "\", which is neither \"1\" nor \"-1\".");

  if(m_writeSystemEvery) { // only then, the following attribute is relevant
    if(m_outputFileName == "undefined")
      throw gError("LinSolveTensorVector::setup", "Attribute 'nameOutputFile' has value \"undefined\".");
  }

  m_colour1 = M_MANAGER->getColour(m_species1);
  m_colour2 = M_MANAGER->getColour(m_species2);

  m_cps[0] = M_MANAGER->cp(m_colour1, m_colour1);
  if(m_colour1 != m_colour2) { // otherwise the entries stay at NULL
    m_cps[1] = M_MANAGER->cp(m_colour1, m_colour2);
    m_cps[2] = M_MANAGER->cp(m_colour2, m_colour2);
  }

  if(Particle::s_tag_format[m_colour1].attrExists(m_solSymbolName)) {
    throw gError("LinSolveTensorVector::setup", "Conflicting symbol \"" + m_solSymbolName + "\" already exists for species " + m_species1 + ".");
  }
  else {
    // Should be OK to be persistent, hence not set to zero. Anyway we WRITE the solution and do not INCREMENT. Persistency should additionally help, that the solution is available as long as no new one is computed.
    m_solSymbolOffset1 = Particle::s_tag_format[m_colour1].addAttribute(m_solSymbolName, DataFormat::POINT, /*persistent=*/true, m_solSymbolName).offset;
  }

  if(m_colour2 != m_colour1) {
    if(Particle::s_tag_format[m_colour2].attrExists(m_solSymbolName)) {
      throw gError("LinSolveTensorVector::setup", "Conflicting symbol \"" + m_solSymbolName + "\" already exists for species " + m_species2 + ".");
    }
    else {
      // Should be OK to be persistent, hence not set to zero. Anyway we WRITE the solution and do not INCREMENT. Persistency should additionally help, that the solution is available as long as no new one is computed.
      m_solSymbolOffset2 = Particle::s_tag_format[m_colour2].addAttribute(m_solSymbolName, DataFormat::POINT, /*persistent=*/true, m_solSymbolName).offset;
    }
  }
  else m_solSymbolOffset2 = m_solSymbolOffset1;

  if(!Particle::s_tag_format[m_colour1].attrExists(m_RHSsymbolName)) {
    throw gError("LinSolveTensorVector::setup", "Symbol \"" + m_RHSsymbolName + "\" not found for species " + m_species1 + ".");
  }
  else {
    m_RHSsymbolOffset1 = Particle::s_tag_format[m_colour1].offsetByName(m_RHSsymbolName);
  }

  if(m_colour2 != m_colour1) {
    if(!Particle::s_tag_format[m_colour2].attrExists(m_RHSsymbolName)) {
      throw gError("LinSolveTensorVector::setup", "Symbol \"" + m_RHSsymbolName + "\" not found for species " + m_species2 + ".");
    }
    else {
      m_RHSsymbolOffset2 = Particle::s_tag_format[m_colour2].offsetByName(m_RHSsymbolName);
    }
  }
  else m_RHSsymbolOffset2 = m_RHSsymbolOffset1;

  int cpIdx = -1;
  for(vector<ColourPair*>::iterator cpit = m_cps.begin(); *cpit != NULL; ++cpit) {
    ++cpIdx;
    ColourPair* cp = *cpit;

    m_functions[cpIdx] = new FunctionPair();
    m_functions[cpIdx]->setExpression(m_expression);
    m_functions[cpIdx]->setColourPair(cp);
    m_functions[cpIdx]->setReturnType(Variant::TENSOR);
    
    m_1stParticleFactors[cpIdx] = new FunctionPair();
    m_1stParticleFactors[cpIdx]->setExpression(m_1stPExpression);
    m_1stParticleFactors[cpIdx]->setColourPair(cp);
    m_1stParticleFactors[cpIdx]->setReturnType(Variant::TENSOR);
    
    m_2ndParticleFactors[cpIdx] = new FunctionPair();
    m_2ndParticleFactors[cpIdx]->setExpression(m_2ndPExpression);
    m_2ndParticleFactors[cpIdx]->setColourPair(cp);
    m_2ndParticleFactors[cpIdx]->setReturnType(Variant::TENSOR);

  }

  M_CONTROLLER->registerForSetupAfterParticleCreation(this);

}

void LinSolveTensorVector::setupAfterParticleCreation() 
{
  size_t nOfParticles = M_PHASE->particles(m_colour1).size();
  if(m_colour1 != m_colour2)
    nOfParticles += M_PHASE->particles(m_colour2).size();

  m_mat.resize(nOfParticles*nOfParticles*SPACE_DIMS_SQUARED);

  // fill in zeroes
  std::fill(m_mat.begin(), m_mat.end(), 0.);
}

void LinSolveTensorVector::call(size_t timestep) {
// MSG_DEBUG("LinSolveTensorVector::call", "START, timestep = " << timestep);
  
  Phase* phase = M_PHASE;
  
  size_t nParticles1 = phase->returnNofPartC(m_colour1);
  size_t nParticles = nParticles1;
//   size_t nParticles2 = nParticles1;
  if(m_colour1 != m_colour2) {
//     nParticles2 = phase->nofFreeParticles(m_colour2);
//     nParticles += nParticles2;
    nParticles += phase->returnNofPartC(m_colour2);
  }
  
  // offsets for accessing m_mat for second and third ColourPair
  // for the three ColourPairs, we will have for (nPOffset1,nPOffset2): (0,0), (0,1), (1,1)
  // see the increments at end of loop
  int nPOffset1 = 0;
  int nPOffset2 = 0;
  int cpIdx = -1;
  for(vector<ColourPair*>::iterator cpit = m_cps.begin(); *cpit != NULL; ++cpit) {
    ++cpIdx;
    FunctionPair* m_function = m_functions[cpIdx];
    FunctionPair* m_1stParticleFactor = m_1stParticleFactors[cpIdx];
    FunctionPair* m_2ndParticleFactor = m_2ndParticleFactors[cpIdx];
    ColourPair* cp = *cpit;
    FOR_EACH_PAIR
      (cp,  
       if(pair->abs() < m_cutoff) {
	 tensor_t temp;
	 
	 // compute the pair-expression
	 (*m_function)(&temp, pair);
	 
	 tensor_t tempFirst;
	 tensor_t tempSecond;
	 // compute the particle-expressions
	 (*m_1stParticleFactor)(&tempFirst, pair);
	 (*m_2ndParticleFactor)(&tempSecond, pair);
	 
	 Particle* first = pair->firstPart();
	 Particle* second = pair->secondPart();
	 // the offsets only become active for the 2nd or latest the 3rd ColourPair 
	 // (see definition above), by assuming that the first rows and columns of 
	 // m_mat are for species1 and the following for species2
	 size_t firstSlot = nPOffset1*nParticles1 + first->mySlot;
	 size_t secondSlot = nPOffset2*nParticles1 + second->mySlot;
	 // MSG_DEBUG("LinSolveTensorVector::compute", "temp = " << temp);
	 
	 // write into matrix
	 // FIXME: think about sparse-matrix format
	 if(pair->actsOnFirst()) {
	   // component-wise product for tensor_t
	   tempFirst *= temp;
	   
	   // precompute
	   size_t firstSlotSPACE_DIMS_SQUARED = firstSlot*SPACE_DIMS_SQUARED;

	   for(size_t i = 0; i < SPACE_DIMS; ++i) {
	     // precompute
	     size_t iSPACE_DIMS = i*SPACE_DIMS;
	     for(size_t j = 0; j < SPACE_DIMS; ++j) {
	       // off-diagonal
	       /* 		   m_mat[j+i*SPACE_DIMS+secondSlot*SPACE_DIMS*SPACE_DIMS+firstSlot*SPACE_DIMS*SPACE_DIMS*nParticles] = tempFirst(i,j); */
	       m_mat[j+iSPACE_DIMS+secondSlot*SPACE_DIMS_SQUARED+firstSlotSPACE_DIMS_SQUARED*nParticles] = tempFirst(i,j);
	       /* 		   m_mat[j+(i+(secondSlot+firstSlot*nParticles)*SPACE_DIMS)*SPACE_DIMS] = tempFirst(i,j); */
	       
	       // subtract from diagonal for first particle
	       /* 		   m_mat[j+i*SPACE_DIMS+firstSlot*SPACE_DIMS*SPACE_DIMS+firstSlot*SPACE_DIMS*SPACE_DIMS*nParticles] -= tempFirst(i,j); */
	       m_mat[j+iSPACE_DIMS+firstSlotSPACE_DIMS_SQUARED*(1+nParticles)] -= tempFirst(i,j);
	       /* 		   m_mat[j+(i+(firstSlot+firstSlot*nParticles)*SPACE_DIMS)*SPACE_DIMS] -= tempFirst(i,j); */
	       /* 		   m_mat[j+(i+(firstSlot*(1+nParticles))*SPACE_DIMS)*SPACE_DIMS] -= tempFirst(i,j); */
	     }
	   }
	   // MSG_DEBUG("LinSolveTensorVector::compute", "AFTER updating first particle entries: NO DEBUG OUTPUT SPECIFIED YET!);
	 }
	 
	 if(pair->actsOnSecond()) {
	   // component-wise product for tensor_t
	   tempSecond *= temp*m_symmetry;
	   
	   // precompute
	   size_t secondSlotSPACE_DIMS_SQUARED = secondSlot*SPACE_DIMS_SQUARED;

	   for(size_t i = 0; i < SPACE_DIMS; ++i) {
	     // precompute
	     size_t iSPACE_DIMS = i*SPACE_DIMS;
	     for(size_t j = 0; j < SPACE_DIMS; ++j) {
	       // i and j can stay the same (REALLY?!?), but first and second swap compared to above
	       // off-diagonal
	       /* 		   m_mat[j+i*SPACE_DIMS+firstSlot*SPACE_DIMS*SPACE_DIMS+secondSlot*SPACE_DIMS*SPACE_DIMS*nParticles] = tempSecond(i,j); */
	       m_mat[j+iSPACE_DIMS+firstSlot*SPACE_DIMS_SQUARED+secondSlotSPACE_DIMS_SQUARED*nParticles] = tempSecond(i,j);
	       /* 		   m_mat[j+(i+(firstSlot+secondSlot*nParticles)*SPACE_DIMS)*SPACE_DIMS] = tempSecond(i,j); */
	       
	       // subtract from diagonal for first particle
	       /* 		   m_mat[j+i*SPACE_DIMS+secondSlot*SPACE_DIMS*SPACE_DIMS+secondSlot*SPACE_DIMS*SPACE_DIMS*nParticles] -= tempSecond(i,j); */
	       m_mat[j+iSPACE_DIMS+secondSlotSPACE_DIMS_SQUARED*(1+nParticles)] -= tempSecond(i,j);
	       /* 		   m_mat[j+(i+(secondSlot+secondSlot*nParticles)*SPACE_DIMS)*SPACE_DIMS] -= tempSecond(i,j); */
	       /* 		   m_mat[j+(i+(secondSlot*(1+nParticles))*SPACE_DIMS)*SPACE_DIMS] -= tempSecond(i,j); */
	     }
	   }	       
	   // MSG_DEBUG("LinSolveTensorVector::compute", "AFTER updating second particle entries: NO DEBUG OUTPUT SPECIFIED YET!);
	 }
       }
       
       ); // end FOR_EACH_PAIR
    
    // for the three ColourPairs, we will have for (nPOffset1,nPOffset2): (0,0), (0,1), (1,1)
    nPOffset1 += nPOffset2;
    nPOffset2 += 1-nPOffset1; 
  } // end of for(vector<ColourPair*>::iterator cpit = m_cps.begin(); *cpit != NULL; ++cpit)
  
  
    //////////////////// write matrix to file //////////////////////
  
  if(m_writeSystemEvery) {
    // is it time to write again?
    if(!m_countdown) {
      // reset countdown
      m_countdown = m_writeSystemEvery-1;
      
      if(m_binary) {
	// open new file for matrix
	m_s.open(make_filename(m_outputFileName, m_step_counter).c_str(), ios::out|ios::binary);
	// write
	m_s.write((char*) &m_mat[0], nParticles*nParticles*SPACE_DIMS_SQUARED*sizeof(double));
	m_s.close();

	// open new file for RHS
	m_s.open(make_filename(append_to_filename(m_outputFileName, "RHS"), m_step_counter).c_str(), ios::out|ios::binary);
	// write
	FOR_EACH_FREE_PARTICLE_C
	  (M_PHASE,m_colour1,
	   m_s.write((char*) &(__iSLFE->tag.pointByOffset(m_RHSsymbolOffset1).coords), SPACE_DIMS*sizeof(double));
	   );
	if(m_colour2 != m_colour1) {
	  FOR_EACH_FREE_PARTICLE_C
	    (M_PHASE,m_colour2,
	     m_s.write((char*) &(__iSLFE->tag.pointByOffset(m_RHSsymbolOffset2).coords), SPACE_DIMS*sizeof(double));
	     );
	}
	m_s.close();

      } // end of if(m_binary)
      else { // so it will be ASCII output
	// open new file for matrix
	m_s.open(make_filename(m_outputFileName, m_step_counter).c_str());
	m_s.flags(ios::scientific);
	writeHeader();
	// The index-order is: p1, p2, i, j (from slowest to fastest index)
	// Hence row p1 in 3D is (p1,0,0,0),(p1,0,0,1),...(p1,0,2,2),(p1,1,0,0),...(p1,nP,2,2),
	// i.e., all entries of matrix for (p1,p2) appear consecutively in ONE row
	// i.e., the ASCII-file has nParticles rows and nParticles*SPACE_DIMS_SQUARED columns
	// One row for each particle
	for(size_t p1 = 0; p1 < nParticles; ++p1) {
	  size_t p1SPACE_DIMS_SQUAREDnP = p1*SPACE_DIMS_SQUARED*nParticles ;
	  for(size_t p2 = 0; p2 < nParticles; ++p2) {
	    size_t p1p2Index = p1SPACE_DIMS_SQUAREDnP + p2*SPACE_DIMS_SQUARED;
	    for(size_t i = 0; i < SPACE_DIMS; ++i) {
	      size_t ip1p2Index = i*SPACE_DIMS+p1p2Index;
	      for(size_t j = 0; j < SPACE_DIMS; ++j) {
		m_s << m_mat[j+ip1p2Index] << " ";
	      }	  
	    }
	  }
	  m_s << endl; // new p1 goes in new row
	}
	m_s.close();
 
	// open new file for RHS
	m_s.open(make_filename(append_to_filename(m_outputFileName, "RHS"), m_step_counter).c_str());
	m_s.flags(ios::scientific);
	FOR_EACH_FREE_PARTICLE_C
	  (M_PHASE,m_colour1,
	   m_s << __iSLFE->tag.pointByOffset(m_RHSsymbolOffset1)<< " ";
	   m_s << endl; // new p goes in new row
	   );
	if(m_colour2 != m_colour1) {
	  FOR_EACH_FREE_PARTICLE_C
	    (M_PHASE,m_colour2,
	     m_s << __iSLFE->tag.pointByOffset(m_RHSsymbolOffset2)<< " ";
	     m_s << endl; // new p goes in new row
	     );
	}
	m_s.close();

     } // end of else of if(m_binary)
      ++m_step_counter;
    }
    else --m_countdown; // count down the countdown
  } // end of if(m_writeSystemEvery)
  
  
  ///////// Solving the linear system ////////////////////////////////////////////
    
    // FIXME: ADD THE SOLVING HERE
    
    
    ///////// Loop again over pairs and particles and set entries back to zero /////
    
    // (i.e., only those that are !=0 instead of brute force setting everything 
    // to zero)
    // loop for off-diagonal entries
    
    // offsets for accessing m_mat for second and third ColourPair
    // for the three ColourPairs, we will have for (nPOffset1,nPOffset2): (0,0), (0,1), (1,1)
    // see the increments at end of loop
    nPOffset1 = 0;
    nPOffset2 = 0;
    for(vector<ColourPair*>::iterator cpit = m_cps.begin(); *cpit != NULL; ++cpit) {
      ColourPair* cp = *cpit;
      FOR_EACH_PAIR
	(cp,  
	 if(pair->abs() < m_cutoff) {	   
	   Particle* first = pair->firstPart();
	   Particle* second = pair->secondPart();
	   // the offsets only become active for the 2nd or latest the 3rd ColourPair 
	   // (see definition above), by assuming that the first rows and columns of 
	   // m_mat are for species1 and the following for species2
	   size_t firstSlot = nPOffset1*nParticles1 + first->mySlot;
	   size_t secondSlot = nPOffset2*nParticles1 + second->mySlot;
	   // MSG_DEBUG("LinSolveTensorVector::compute", "temp = " << temp);
	   
	   // write into matrix
	   // FIXME: think about sparse-matrix format
	   if(pair->actsOnFirst()) {
	     // precompute
	     size_t firstSlotSPACE_DIMS_SQUARED = firstSlot*SPACE_DIMS_SQUARED;
	     
	     for(size_t i = 0; i < SPACE_DIMS; ++i) {
	       // precompute
	       size_t iSPACE_DIMS = i*SPACE_DIMS;
	       for(size_t j = 0; j < SPACE_DIMS; ++j) {
		 // off-diagonal
		 /* 		   m_mat[j+i*SPACE_DIMS+secondSlot*SPACE_DIMS*SPACE_DIMS+firstSlot*SPACE_DIMS*SPACE_DIMS*nParticles] = 0; */
		 m_mat[j+iSPACE_DIMS+secondSlot*SPACE_DIMS_SQUARED+firstSlotSPACE_DIMS_SQUARED*nParticles] = 0;
		 /* 		   m_mat[j+(i+(secondSlot+firstSlot*nParticles)*SPACE_DIMS)*SPACE_DIMS] = 0; */
		 
	       }
	     }
	     // MSG_DEBUG("LinSolveTensorVector::compute", "AFTER updating first particle entries: NO DEBUG OUTPUT SPECIFIED YET!);
	   }
	   
	   if(pair->actsOnSecond()) {
	     // precompute
	     size_t secondSlotSPACE_DIMS_SQUARED = secondSlot*SPACE_DIMS_SQUARED;
	     
	     for(size_t i = 0; i < SPACE_DIMS; ++i) {
	       // precompute
	       size_t iSPACE_DIMS = i*SPACE_DIMS;
	       for(size_t j = 0; j < SPACE_DIMS; ++j) {
		 // i and j can stay the same (REALLY?!?), but first and second swap compared to above
		 // off-diagonal
		 /* 		   m_mat[j+i*SPACE_DIMS+firstSlot*SPACE_DIMS*SPACE_DIMS+secondSlot*SPACE_DIMS*SPACE_DIMS*nParticles] = 0; */
		 m_mat[j+iSPACE_DIMS+firstSlot*SPACE_DIMS_SQUARED+secondSlotSPACE_DIMS_SQUARED*nParticles] = 0;
		 /* 		   m_mat[j+(i+(firstSlot+secondSlot*nParticles)*SPACE_DIMS)*SPACE_DIMS] = 0; */
		 
	       }
	     }	       
	     // MSG_DEBUG("LinSolveTensorVector::compute", "AFTER updating second particle entries: NO DEBUG OUTPUT SPECIFIED YET!);
	   }
	 }       
	 ); // end FOR_EACH_PAIR
      
      // for the three ColourPairs, we will have for (nPOffset1,nPOffset2): (0,0), (0,1), (1,1)
      nPOffset1 += nPOffset2;
      nPOffset2 += 1-nPOffset1; 

    } // end of for(vector<ColourPair*>::iterator cpit = m_cps.begin(); *cpit != NULL; ++cpit)
    
    
    // loop for diagonal entries, first species
    FOR_EACH_PARTICLE_C
      (phase, m_colour1,
       size_t slot = __iSLFE->mySlot;
       // precompute
       size_t slotSPACE_DIMS_SQUARED = slot*SPACE_DIMS_SQUARED;
       for(size_t i = 0; i < SPACE_DIMS; ++i) {
	 // precompute
	 size_t iSPACE_DIMS = i*SPACE_DIMS;
	 for(size_t j = 0; j < SPACE_DIMS; ++j) {	   
	   /* 	   m_mat[j+i*SPACE_DIMS+slot*SPACE_DIMS*SPACE_DIMS+slot*SPACE_DIMS*SPACE_DIMS*nParticles] = 0; */
	   m_mat[j+iSPACE_DIMS+slotSPACE_DIMS_SQUARED*(1+nParticles)] = 0;
	   /* 	   m_mat[j+(i+(slot+slot*nParticles)*SPACE_DIMS)*SPACE_DIMS] = 0; */
	   /* 	   m_mat[j+(i+(slot*(1+nParticles))*SPACE_DIMS)*SPACE_DIMS] = 0; */
	 }
       }
       );
    
    // loop for diagonal entries, second species
    if(m_colour1 != m_colour2) {
      FOR_EACH_PARTICLE_C
	(phase, m_colour2,
	 size_t slot = __iSLFE->mySlot;
	 // precompute
	 size_t slotSPACE_DIMS_SQUARED = slot*SPACE_DIMS_SQUARED;
	 for(size_t i = 0; i < SPACE_DIMS; ++i) {
	   // precompute
	   size_t iSPACE_DIMS = i*SPACE_DIMS;
	   for(size_t j = 0; j < SPACE_DIMS; ++j) {	   
	     /* 	   m_mat[j+i*SPACE_DIMS+slot*SPACE_DIMS*SPACE_DIMS+slot*SPACE_DIMS*SPACE_DIMS*nParticles] = 0; */
	     m_mat[j+iSPACE_DIMS+slotSPACE_DIMS_SQUARED*(1+nParticles)] = 0;
	     /* 	   m_mat[j+(i+(slot+slot*nParticles)*SPACE_DIMS)*SPACE_DIMS] = 0; */
	     /* 	   m_mat[j+(i+(slot*(1+nParticles))*SPACE_DIMS)*SPACE_DIMS] = 0; */
	   }
	 }
	 );
    }
}

void LinSolveTensorVector::writeHeader() {
  m_s << 
    "#Matrix written out by LinSolveTensorVector.\n"
     "#The index-order is: p1, p2, i, j (from slowest to fastest index)\n"
     "#Hence row p1 in 3D is (p1,0,0,0),(p1,0,0,1),...(p1,0,2,2),(p1,1,0,0),...(p1,nP,2,2),\n"
     "#i.e., all entries of matrix for (p1,p2) appear consecutively in ONE row\n"
    "#i.e., the ASCII-file has nParticles rows and nParticles*SPACE_DIMS_SQUARED columns" << endl;
}
