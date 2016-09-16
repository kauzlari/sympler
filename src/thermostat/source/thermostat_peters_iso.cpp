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



#include "thermostat_peters_iso.h"

#include "phase.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"
#include "random.h"
#include "pair_creator.h"
// #include "node.h"

using namespace std;


#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()

const Callable_Register<ThermostatPetersIso> thermostat_peters_iso("ThermostatPetersIso");


ThermostatPetersIso::ThermostatPetersIso(Simulation* sim)
  : ThermostatPeters(sim)
{
  init();
}


void ThermostatPetersIso::init()
{
  m_properties.setClassName("ThermostatPetersIso");

  m_properties.setDescription(
    "Isothermal Peters thermostat. See: E. A. J. F. Peters, Europhys. Lett. 66 (3), pp. 311.\n"
    "NOTE: Equal masses m=1 are assumed for all thermalised particles!"

);

  DOUBLEPC
    (kBToverM, m_temperature, 0,
     "k_BT to thermalize to.");

  DOUBLEPC
    (mass1, m_mass1, 0,
     "Particle mass of species 1.");
  DOUBLEPC
    (mass2, m_mass2, 0,
     "Particle mass of species 2.");

  m_temperature = 1;
  m_mass1 = 1;
  m_mass2 = 1; 
}


void ThermostatPetersIso::setup()
{
  ThermostatPeters::setup();
  m_massred = m_mass1*m_mass2/(m_mass1+m_mass2);
}


void ThermostatPetersIso::thermalize(Phase* phase)
{
  phase->pairCreator()->createDistances();

  FOR_EACH_PAIR__PARALLEL
      (ThermostatPetersIso,
       m_cp,
       if(pair->abs() < self->m_cutoff) {
       
         double temp;
         self->m_dissipation(&temp, &(*pair));
         point_t g;
         g.assign(0);
         double weight = self->m_wf->interpolate(pair, g);
//       RandomNumberGenerator m_rng;
         /* ---- Peters Scheme II ---- */

         weight *= weight * temp/*self->m_dissipation*/ * self->m_dt / m_massred;

         double a = (1-exp(-weight)) * m_massred;
         double b = sqrt
               ((1-exp(-2*weight))*self->m_temperature *  m_massred)*m_rng.normal(1);


         /* ---- Peter Scheme I ---- */
       /*
         double a = m_dissipation*weight*weight*m_dt;
         double b = sqrt(2*m_temperature*
         a*(1-a))*m_rng.normal(1);
       */

#ifdef ENABLE_PTHREADS
       pair->firstPart()->lock();
       pair->secondPart()->lock();
#endif

       point_t vij = pair->firstPart()->v - pair->secondPart()->v;
       point_t rij = pair->cartesian()/pair->cartesian().abs();

       double vr = vij*rij;
       double h = -a*vr + b;

       point_t dv = h*rij;

       if (pair->actsOnFirst()) {
         pair->firstPart()->v += dv/m_mass1;
       }

       if (pair->actsOnSecond()) {
         pair->secondPart()->v -= dv/m_mass2;
       }

#ifdef ENABLE_PTHREADS
       pair->secondPart()->unlock();
       pair->firstPart()->unlock();
#endif
       }
      );
}
