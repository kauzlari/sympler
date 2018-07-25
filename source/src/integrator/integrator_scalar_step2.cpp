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


#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "particle.h"
#include "controller.h"
#include "simulation.h"
#include "integrator_scalar_step2.h"

using namespace std;


#define M_CONTROLLER  ((Controller*) m_parent)
#define M_SIMULATION  ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE  M_SIMULATION->phase()


const Integrator_Register<IntegratorScalarStep2> integrator_scalar_step2("IntegratorScalarStep2");

//---- Constructors/Destructor ----

IntegratorScalarStep2::IntegratorScalarStep2(Controller *controller): IntegratorScalar(controller)
{
  init();
}


IntegratorScalarStep2::~IntegratorScalarStep2()
{
}

//---- Methods ----

void IntegratorScalarStep2::init()
{
  m_properties.setClassName("IntegratorScalarStep2");

  m_properties.setDescription
    (
     "Adds an additional scalar degree of freedom, to the particles of "
     "the specified colour and integrates it according to the following "
     "scheme:\n"
     "\n"
     "integration-step1: no activity\n"
     "integration-step2: s(t + dt) = s(t) + dt * F(t)\n"
     "\n"
     "Here, s is the scalar, F is its flux (usually computed with Force "
     "modules such as FPairScalar or FParticleScalar), t is time and dt is "
     "the size of the integration time step (defined in the Controller). "
     "Further information on the integration-steps, including their place "
     "in the total SYMPLER workflow, can be found with the help option "
     "\"--help workflow\"."
     );

}

void IntegratorScalarStep2::setup()
{
  IntegratorScalar::setup();
}


void IntegratorScalarStep2::integrateStep1()
{
}


void IntegratorScalarStep2::integrateStep2()
{
  Phase *phase = M_PHASE;

  size_t force_index = M_CONTROLLER->forceIndex();

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, m_colour, this,
     
     i->tag.doubleByOffset(((IntegratorScalarStep2*) data)->m_scalar_offset)
     +=
     ((IntegratorScalarStep2*) data)->m_dt *
     (
      i->tag.doubleByOffset(((IntegratorScalarStep2*) data)
			    -> m_force_offset[force_index]
			    )
      );
     );
}
