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


#include "integrator_ISPH_const_rho_FDRHS.h"

#include "gen_f.h"
#include "phase.h"
#include "threads.h"
#include "controller.h"
#include "simulation.h"
#include "cell.h"

using namespace std;


#define M_CONTROLLER ((Controller*) m_parent)
#define M_SIMULATION ((Simulation*) M_CONTROLLER->parent())
#define M_PHASE M_SIMULATION->phase()

#define M_MANAGER M_PHASE->manager()
const Integrator_Register<IntegratorISPHconstRhoFDRHS> integrator_ISPH_const_rho_FDRHS("IntegratorISPHconstRhoFDRHS");


//---- Constructors/Destructor ----

IntegratorISPHconstRhoFDRHS::IntegratorISPHconstRhoFDRHS(Controller *controller):IntegratorISPHconstRho(controller)
{
  init();
}


IntegratorISPHconstRhoFDRHS::~IntegratorISPHconstRhoFDRHS()
{
}


//---- Methods ----

void IntegratorISPHconstRhoFDRHS::init()
{
  // some modules need to know whether there is an Integrator,
  // which changes positions, that's why the following
  m_properties.setClassName("IntegratorPosition");
  m_properties.setName("IntegratorISPHconstRhoFDRHS");

  m_properties.setDescription
    (m_properties.description() +
     "\nDefinition of the right-hand-side (RHS) of the PPE for "
     "IntegratorISPHconstRhoFDRHS:\n"
     "RHS = (rho_0-rho_adv)/rho_0/dt^2\n"
     "where rho_0 is the constant reference density, rho_adv is the "
     "advected density resulting from all forces except the pressure "
     "gradient forces, and dt is the integration time step. "
     );
  
}

void IntegratorISPHconstRhoFDRHS::addRHStoNewPressure(size_t colour) {

  Phase* phase = M_PHASE;
  double rhsDenom = 1./(m_dt*m_dt*m_rho0);
  size_t advDensityOffset = m_advDensityOffset[colour];
  size_t pressureIterNewOffset = m_pressureIterNewOffset[colour];

  FOR_EACH_FREE_PARTICLE_C__PARALLEL
    (phase, colour, this,
     const Data& pTag = i->tag;
     double& newP = pTag.doubleByOffset(pressureIterNewOffset);
     // add RHS  
     newP += rhsDenom*(m_rho0 - pTag.doubleByOffset(advDensityOffset));

     MSG_DEBUG("IntegratorISPHconstRhoFDRHS::addRHStoNewPressure", "p=" << i->mySlot << ", RHS=" << rhsDenom*(m_rho0 - pTag.doubleByOffset(advDensityOffset)));

     );

}
