/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "particle.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_vector_lambda.h"

using namespace std;


#define M_CONTROLLER  ((Controller*) m_parent)
#define M_SIMULATION  ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE  M_SIMULATION->phase()


const Integrator_Register<IntegratorVectorLambda> integrator_vector_lambda("IntegratorVectorLambda");

//---- Constructors/Destructor ----

IntegratorVectorLambda::IntegratorVectorLambda(Controller *controller): IntegratorVector(controller)
{
  init();
}


IntegratorVectorLambda::~IntegratorVectorLambda()
{
}

//---- Methods ----

void IntegratorVectorLambda::init()
{

  m_properties.setClassName("IntegratorVectorLambda");

  m_properties.setDescription
    (
     "Adds an additional vector degree of freedom, "
     "to the particles specified and integrates it with a second "
     "order accurate predictor-corrector scheme."
     "\nPredictor step: predicted = old + lambda*dt*f_new"
     "\nCorrector step: new = predicted + 0.5*dt*(f_new - f_old)"
     "\nUsually, the forces are updated between the two steps. "
     );

  DOUBLEPC
    (lambda, m_lambda, 0,
     "Lambda parameter for adjustement of predictor-corrector scheme.");

  m_lambda = 0.5;
}

void IntegratorVectorLambda::setup()
{
  IntegratorVector::setup();
}


void IntegratorVectorLambda::integrateStep1()
{
  Phase *phase = M_PHASE;

  size_t force_index = M_CONTROLLER->forceIndex();

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     
     for(size_t j = 0; j < SPACE_DIMS; ++j) {  
       
       i->tag.pointByOffset(((IntegratorVectorLambda*) data)
			    -> m_vector_offset)[j] +=
	 
	 ((IntegratorVectorLambda*) data)->m_dt *
	 
	 m_lambda *
	 i->tag.pointByOffset(((IntegratorVectorLambda*) data) 
			      -> m_force_offset[force_index])[j];
     }
     );
}


void IntegratorVectorLambda::integrateStep2()
{

  Phase *phase = M_PHASE;

  size_t force_index = M_CONTROLLER->forceIndex();
  size_t other_force_index = (force_index+1)&(FORCE_HIST_SIZE-1);

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     
     i->tag.pointByOffset(((IntegratorVectorLambda*) data)
			  ->m_vector_offset) +=
     ((IntegratorVectorLambda*) data)->m_dt *
     0.5*(
	  // new flux
	  i->tag.pointByOffset(((IntegratorVectorLambda*) data)
				 -> m_force_offset[force_index]
				 )
	  -
	  // old flux
	  i->tag.pointByOffset(((IntegratorVectorLambda*) data)
				 -> m_force_offset[other_force_index]
				 )    
	  );
     );  
  
}
