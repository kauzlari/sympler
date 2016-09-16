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



#include "thermostat_peters_ec.h"


#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"
#include "random.h"
#include "pair_creator.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const Callable_Register<ThermostatPetersEnergyConserving> thermostat_peters_ec("ThermostatPetersEnergyConserving");


ThermostatPetersEnergyConserving::ThermostatPetersEnergyConserving(Simulation* sim)
  : ThermostatPeters(sim)
{
  init();
}


void ThermostatPetersEnergyConserving::init()
{
  m_properties.setClassName("ThermostatPetersEnergyConserving");

  m_properties.setDescription(
    "Modified Peters thermostat including energy conservation. \n"
    "See: Pastewka, L.; Kauzlaric, D.; Greiner, A. & Korvink, J. G.; Phys. Rev. E, 2006, 73, 037701\n"
    "NOTE: Equal masses m=1 are assumed for all thermalised particles!"
  );
}


void ThermostatPetersEnergyConserving::setup()
{
  ThermostatPeters::setup();

  m_ie.first = 
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species.first);
  m_ie.second = 
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species.second);

  if (!m_ie.first && !m_ie.second)
    throw gError
      ("TheromstatPetersEnergyConserving::read",
       "You cannot use this object without IntegratorEnergy for the"
           " corresponding species.");

  m_energy_offset.first = Particle::s_tag_format[m_cp->firstColour()].attrByName("internal_energy").offset;
  m_energy_offset.second = Particle::s_tag_format[m_cp->secondColour()].attrByName("internal_energy").offset;
}


void ThermostatPetersEnergyConserving::thermalize(Phase* phase)
{
  phase->pairCreator()->createDistances();

//     // particle tracking
//     FOR_EACH_FREE_PARTICLE_C
//         (M_SIMULATION->phase(), 0,
//          if(i->mySlot == 183)
//          {
//            MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", i->mySlot << " before(N=" << M_SIMULATION->phase()->returnNofPart() << "):"
//                << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1]);
//          }
//         );

  
  /* Loop over all pairs and apply the update rule. Note that only
     the interpolate method of the weighting function is used.
  */

  FOR_EACH_PAIR__PARALLEL
    (ThermostatPetersEnergyConserving,
     m_cp,
     if(pair->abs() < self->m_cutoff) {
       
       
/*       // particle tracking
       if(pair->firstPart()->mySlot == 183 && pair->firstPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", "183: v_before = " << pair->firstPart()->v << endl << "ie_before = " << pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) << endl << "T = " << pair->firstPart()->tag.doubleByOffset(24) << endl << "rho_lucy1_before = " << pair->firstPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_before = " << pair->firstPart()->tag.tensorByOffset(40) << endl << "PARTNER: c = " << pair->secondPart()->c << ", v_before = " << pair->secondPart()->v << endl << "ie_before = " << pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) << endl << "rho_lucy1_before = " << pair->secondPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_before = " << pair->secondPart()->tag.tensorByOffset(40));
       if(pair->secondPart()->mySlot == 183 && pair->secondPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", "183: v_before = " << pair->secondPart()->v << endl << "ie_before = " << pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) << endl << "T = " << pair->secondPart()->tag.doubleByOffset(24) << endl << "rho_lucy1_before = " << pair->secondPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_before = " << pair->secondPart()->tag.tensorByOffset(40) << endl << "PARTNER: c = " << pair->firstPart()->c << ", v_before = " << pair->firstPart()->v << endl << "ie_before = " << pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) << endl << "rho_lucy1_before = " << pair->firstPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_before = " << pair->firstPart()->tag.tensorByOffset(40));*/
       
       
       double temp;
       self->m_dissipation(&temp, &(*pair));
       point_t g;
       g.assign(0);
       double weight = self->m_wf->interpolate(pair, g);
       //RandomNumberGenerator m_rng;

       
//        // particle tracking
//        if(pair->firstPart()->mySlot == 183 && pair->firstPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", "183: temp = " << temp << ", weight = " << weight);
//        if(pair->secondPart()->mySlot == 183 && pair->secondPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize",  "183: temp = " << temp << ", weight = " << weight);

       
              
       /* ---- Peters Scheme II ---- */
       weight *= 2 * temp/*self->m_dissipation*/ * self->m_dt;

#ifdef ENABLE_PTHREADS
       pair->firstPart()->lock();
       pair->secondPart()->lock();
#endif

       /*-------*/

       double a = (1-exp(-weight))/2;
       double b = m_rng.normal
         (sqrt(
	       (1-exp(-2*weight))
	       /
	       (self->m_ie.first->reciprocalTemperature(*pair->firstPart())
		+
		self->m_ie.second->reciprocalTemperature(*pair->secondPart()))
	       ));
    
       
       
//        if(pair->firstPart()->mySlot == 183 && pair->firstPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", "183: a = " << a << ", b = " << b);
//        if(pair->secondPart()->mySlot == 183 && pair->secondPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize",  "183: a = " << a << ", b = " << b);

       
       
       /* ---- Peter Scheme I ----

       double a = m_dissipation*weight*weight*m_dt;
       double b = sqrt(
       4/(m_ie->reciprocalTemperature(*pair->firstPart())+
       m_ie->reciprocalTemperature(*pair->secondPart()))*
       a*(1-a*m_dt))*m_rng.normal(1);

       */

       point_t vij = pair->firstPart()->v - pair->secondPart()->v;
       point_t eij = pair->cartesian()/pair->cartesian().abs();

       double ve = vij*eij;
       double h = -a*ve + b;

       point_t dv = h*eij;

       double de = h*(h*eij*eij+ve)/2;

       if (pair->actsOnFirst()) {
         pair->firstPart()->v += dv;
         pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) -= de;
       }

       if (pair->actsOnSecond()) {
         pair->secondPart()->v -= dv;
         pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) -= de;
       }

#ifdef ENABLE_PTHREADS
       pair->secondPart()->unlock();
       pair->firstPart()->unlock();
#endif
     }
     
//      // particle tracking
//      if(pair->firstPart()->mySlot == 183 && pair->firstPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", "183: v_after = " << pair->firstPart()->v << endl << "ie_after = " << pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) << endl << "rho_lucy1_after = " << pair->firstPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_after = " << pair->firstPart()->tag.tensorByOffset(40) << endl << "PARTNER: c = " << pair->secondPart()->c << ", v_after = " << pair->secondPart()->v << endl << "ie_after = " << pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) << endl << "rho_lucy1_after = " << pair->secondPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_after = " << pair->secondPart()->tag.tensorByOffset(40));
//      if(pair->secondPart()->mySlot == 183 && pair->secondPart()-> c == 0) MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", "183: v_after = " << pair->secondPart()->v << endl << "ie_after = " << pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) << endl << "rho_lucy1_after = " << pair->secondPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_after = " << pair->secondPart()->tag.tensorByOffset(40) << endl << "PARTNER: c = " << pair->firstPart()->c << ", v_after = " << pair->firstPart()->v << endl << "ie_after = " << pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) << endl << "rho_lucy1_after = " << pair->firstPart()->tag.doubleByOffset(32) << endl << "shear_lucy1_after = " << pair->firstPart()->tag.tensorByOffset(40));

    
    );

//     // particle tracking
//      FOR_EACH_FREE_PARTICLE_C
//          (M_SIMULATION->phase(), 0,
//           if(i->mySlot == 183)
//           {
//             MSG_DEBUG("ThermostatPetersEnergyConserving::thermalize", i->mySlot << " END(N=" << M_SIMULATION->phase()->returnNofPart() << "):"
//                 << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1]);
//           }
//          );

}

