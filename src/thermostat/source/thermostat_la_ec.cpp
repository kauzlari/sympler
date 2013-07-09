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



#include "thermostat_la_ec.h"


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

const Callable_Register<ThermostatLAEnergyConserving> thermostat_la_ec("ThermostatLAEnergyConserving");


ThermostatLAEnergyConserving::ThermostatLAEnergyConserving(Simulation* sim)
  : ThermostatWithRng(sim)
{
  init();
}


void ThermostatLAEnergyConserving::init()
{
  m_properties.setClassName("ThermostatLAEnergyConserving");

  m_properties.setDescription(
    "Lowe-Andersen thermostat. See: C._P._Lowe, Europhys. Lett. 47 (2), pp. 145. Modified to "
    "incorporate energy conservation."
  );

//   double dt;
//   if (M_SIMULATION)
//     dt = M_CONTROLLER->dt();
//   else
//     dt = 0.04;

  STRINGPC
    (species1, m_species.first,
     "First species this thermostat should work on.");

  STRINGPC
    (species2, m_species.second,
     "Second species this thermostat should work on.");

//   m_properties.addProperty
//     ("dissipation", PropertyList::DOUBLE, &m_dissipation,
//      new PLCDoubleRange(0, 1/dt),
//      "Dissipation constant. Propability of pair-thermalization is given by "
//      "dissipation*dt.");

  FUNCTIONPAIRPC
      (probability, m_probability/*, 0*/,
      "Probability of thermalising a pair. \nYou must assure that it always lies in the range [0,1]! \nType some nonsense to "
          "obtain an incomplete list of possible variables and constants. "
          "The expression may contain vectors and tensors,"
          " but as a whole it must represent a scalar.");

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Only pairs within this cut-off radius are considered for thermalization.");

  m_species.first = "UNDEF";
  m_species.second = "UNDEF";
//   m_dissipation = 0.5;
//   m_probability = m_dissipation * dt;
  m_cutoff = -1;
}


void ThermostatLAEnergyConserving::setup()
{
  ThermostatWithRng::setup();

  /* following factor is in general = sqrt(2*k*T/m); currently k=m=1 and m=const. */
//   m_probability = m_dissipation * M_CONTROLLER->dt();

  if(m_cutoff == -1)
    throw gError("ThermostatLAEnergyConserving::read: Cutoff not defined.");

  m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

  m_cp->setCutoff(m_cutoff);
  m_cp->setNeedPairs(true);

  m_ie.first = 
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species.first);
  m_ie.second = 
    (IntegratorEnergy*) M_SIMULATION->controller()->findIntegrator("IntegratorEnergy", m_species.second);

  if (!m_ie.first || !m_ie.second)
    throw gError
      ("TheromstatLAEnergyConserving::read",
       "You cannot use this object without IntegratorEnergy for the"
           " corresponding species.");

  m_energy_offset.first = Particle::s_tag_format[m_cp->firstColour()].attrByName("internal_energy").offset;
  m_energy_offset.second = Particle::s_tag_format[m_cp->secondColour()].attrByName("internal_energy").offset;

  m_probability.setReturnType(Variant::SCALAR);
  m_probability.setColourPair(m_cp);

}


void ThermostatLAEnergyConserving::thermalize(Phase* phase)
{
  phase->pairCreator()->createDistances();
//  RandomNumberGenerator m_rng;
  
  FOR_EACH_PAIR__PARALLEL
    (ThermostatLAEnergyConserving,
     m_cp,

     if(pair -> abs() < self->m_cutoff) {
       /* each pair is rescaled with probability 'm_probability' */
       double prob;
       self->m_probability(&prob, &(*pair));
       if(m_rng.uniform() < prob/*self->m_probability*/)
       {
         double maxwell_factor;
         double newAbsRelVel;
         double energy1;
         double energy2;
         double delta_e;
         point_t delta;

         /* First, we need to compute the maxwell factor for this pair.
            The question remains how to compute the average T
                this choice: Tavg = 2/(1/T1+1/T2)
            like in DPDE.
            Fixme!!! Find proof or evidence. */
/*         maxwell_factor =
           sqrt(4/(m_ie->reciprocalTemperature(*pair->firstPart())+
           m_ie->reciprocalTemperature(*pair->secondPart())));*/

#ifdef ENABLE_PTHREADS
         pair->firstPart()->lock();
         pair->secondPart()->lock();
#endif

//         maxwell_factor = 
//           sqrt(self->m_ie.first->temperature(*pair->firstPart())+self->m_ie.second->temperature(*pair->secondPart()));

         maxwell_factor = sqrt
           (4/
            (self->m_ie.first->reciprocalTemperature(*pair->firstPart())
             +
             self->m_ie.second->reciprocalTemperature(*pair->secondPart()))
           );

         /* newAbsRelVel is the absolute velocity parallel to the line of centres. */
         newAbsRelVel = m_rng.normal(1) * maxwell_factor * pair->abs();

         newAbsRelVel -= (pair->firstPart()->v - pair->secondPart()->v) * pair->cartesian();

         /* after next we have the scalar 
            (0.5/pair->abs())*((v_new_ij dot e_ij) - (v_ij dot e_ij))
            i.e., only multiplication with vector r_ij is missing */
         newAbsRelVel *= 0.5 / pair -> absSquare();

         /* now, delta[i] is computed and added to the old velocities */
         delta = newAbsRelVel * pair->cartesian();

         energy1 = 0.5*(pair->firstPart()->v.absSquare()+pair->secondPart()->v.absSquare());

         /* it acts at least on one, right? */
         if (pair->actsOnFirst())
           pair->firstPart()->v += delta;

         if (pair->actsOnSecond())
           pair->secondPart()->v -= delta;

         energy2 = 0.5*(pair->firstPart()->v.absSquare()+pair->secondPart()->v.absSquare());

         delta_e = (energy2-energy1)/2;

         /* Now, adjust the energy. */
         if (pair->actsOnFirst()) {
           pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) -= delta_e;
           if (pair->firstPart()->tag.doubleByOffset(self->m_energy_offset.first) <= 0)
             throw gError("ThermostatLAEnergyConserving::thermalize", "Energy <= 0. Choose a different s(e)!");
         }

         if (pair->actsOnSecond()) {
           pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) -= delta_e;
           if (pair->secondPart()->tag.doubleByOffset(self->m_energy_offset.second) <= 0)
             throw gError("ThermostatLAEnergyConserving::thermalize", "Energy <= 0. Choose a different s(e)!");
         }
       }

#ifdef ENABLE_PTHREADS
       pair->secondPart()->unlock();
       pair->firstPart()->unlock();
#endif
     }	
  );
}
