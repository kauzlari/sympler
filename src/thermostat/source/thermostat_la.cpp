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




#include "thermostat_with_rng.h"

#include "phase.h"
#include "random.h"
#include "threads.h"
#include "simulation.h"
#include "manager_cell.h"
#include "thermostat_la.h"
#include "pair_creator.h"

using namespace std;

#define M_SIMULATION ((Simulation*) m_parent)
#define M_CONTROLLER M_SIMULATION->controller()
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_PAIRCREATOR M_PHASE->pairCreator()


const Callable_Register<ThermostatLA> thermostat_la("ThermostatLA");



/* --- ThermostatLA --- */

ThermostatLA::ThermostatLA(Simulation* sim)
  : ThermostatWithRng(sim)
{
  init();
}


ThermostatLA::~ThermostatLA()
{
}


void ThermostatLA::init()
{
  m_properties.setClassName("ThermostatLA");

  m_properties.setDescription
    ("Lowe-Andersen thermostat. See: C._P._Lowe, Europhys. Lett. 47 (2), pp. 145. The complete expression added to particle i and subtracted from particle j is: Deltaij'=particleFactor*(Deltaij+particleAddend), where Deltaij is the original term from Lowe's paper.");

  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Only pairs within this cut-off radius are considered for thermalization.");

  m_cutoff = 1;

//   m_properties.addProperty
//     ("probability", PropertyList::DOUBLE, &m_probability,
//      new PLCDoubleRange(0, 1.), "Probability of pair thermalisation");

//   m_probability = 0.5;

  FUNCTIONPAIRPC
      (probability, m_probability/*, 0*/,
      "Probability of thermalising a pair. \nYou must assure that it always lies in the range [0,1]! \nType some nonsense to "
          "obtain an incomplete list of possible variables and constants. "
          "The expression may contain vectors and tensors,"
          " but as a whole it must represent a scalar.");

  FUNCTIONPAIRPC
      (particleFactor_i, m_firstFac/*, 0*/,
      "First particle factor multiplying (piecewise) the relative velocity drawn by the thermostat and the particleAddend. \nType some nonsense to "
          "obtain an incomplete list of possible variables and constants. "
          "The expression may contain vectors and tensors,"
          " but as a whole it must represent a vector.");

  FUNCTIONPAIRPC
      (particleFactor_j, m_secondFac/*, 0*/,
      "Second particle factor multiplying (piecewise) the relative velocity drawn by the thermostat and the particleAddend. \nType some nonsense to "
          "obtain an incomplete list of possible variables and constants. "
          "The expression may contain vectors and tensors,"
          " but as a whole it must represent a vector.");

  FUNCTIONPAIRPC
      (particleAddend_i, m_firstAdd/*, 0*/,
      "First particle addend that is added (piecewise) to the relative velocity drawn by the thermostat. \nType some nonsense to "
          "obtain an incomplete list of possible variables and constants. "
          "The expression may contain vectors and tensors,"
          " but as a whole it must represent a vector.");

  FUNCTIONPAIRPC
      (particleAddend_j, m_secondAdd/*, 0*/,
      "Second particle addend that is added (piecewise) to the relative velocity drawn by the thermostat . \nType some nonsense to "
          "obtain an incomplete list of possible variables and constants. "
          "The expression may contain vectors and tensors,"
          " but as a whole it must represent a vector.");


  STRINGPC
    (species1, m_species.first, "First species of the pair, the Thermostat should thermalise.");
  STRINGPC
    (species2, m_species.second, "Second species of the pair, the Thermostat should thermalise.");

  m_species.first = "UNDEF";
  m_species.second = "UNDEF";
  
  DOUBLEPC
    (temperature, m_temperature, 0,
     "Temperature to thermalize to.");
  m_temperature = 1;
}


void ThermostatLA::setup()
{
  ThermostatWithRng::setup();

  if(m_species.first == "UNDEF") 
    throw gError("ThermostatLA::CPData::read: 'species1' undefined");
  if(m_species.second == "UNDEF") 
    throw gError("ThermostatLA::CPData::read: 'species2' undefined");
  if(m_cutoff < 0) throw gError("ThermostatLA::CPData::read: no positive cutoff defined");
  MSG_DEBUG("ThermostatLA::CPData::read", "m_cutoff=" << m_cutoff);

//   m_cp = ((Simulation*) m_parent)->phase()->manager()->cp(m_species);

  m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);

    
  m_cp->setCutoff(m_cutoff);

  m_cp->setNeedPairs(true);
			
  // following factor is in general = sqrt(2*k*T/m); currently k=m=1 and m=const.
  m_MaxwellFactor = sqrt(2*m_temperature);

  ((Simulation*) m_parent)->maxCutoff = 
    max(((Simulation*) m_parent)->maxCutoff, m_cutoff);

  m_probability.setReturnType(Variant::SCALAR);
  m_probability.setColourPair(m_cp);

  m_firstFac.setReturnType(Variant::VECTOR);
  m_firstFac.setColourPair(m_cp);

  m_secondFac.setReturnType(Variant::VECTOR);
  m_secondFac.setColourPair(m_cp);

  m_firstAdd.setReturnType(Variant::VECTOR);
  m_firstAdd.setColourPair(m_cp);

  m_secondAdd.setReturnType(Variant::VECTOR);
  m_secondAdd.setColourPair(m_cp);

}


void ThermostatLA::thermalize(Phase* phase)
{
  // next line is currently (02/20/2004) not necessary but maybe later
  phase -> pairCreator() -> createDistances();

  FOR_EACH_PAIR__PARALLEL
    (ThermostatLA,
     m_cp, 

     if(pair -> abs() < self->m_cutoff)
       {	
         double prob;
         self->m_probability(&prob, &(*pair));
	 if(m_rng.uniform() < prob/*self->m_probability*/)
	   {
	     double newAbsRelVel = m_rng.normal(1) * self->m_MaxwellFactor * pair->abs();

/*
 MSG_DEBUG("ThermostatLA::thermalize", "BEFORE  pair distance = " << pair->abs() << endl << "velocity  1st PART!!! = " << pair -> firstPart() -> v << endl << "slot  1st PART!!! = " << pair -> firstPart() -> mySlot << endl << "velocity  2nd PART!!! = " << pair -> secondPart() -> v << endl << "slot  2nd PART!!! = " << pair -> secondPart() -> mySlot);*/


#ifdef ENABLE_PTHREADS
	     pair->firstPart()->lock();
	     pair->secondPart()->lock();
#endif
	     // newAbsRelVel = ((v_new_ij dot r_ij) - (v_ij dot r_ij)) 
	     for(size_t i = 0; i < SPACE_DIMS; ++i)
	       {
		 newAbsRelVel -= (pair -> firstPart() -> v[i] - pair -> secondPart() -> 
				  v[i]) * (*pair)[i];
	       }
	     // after next we have the scalar 
	     // (0.5/pair->abs())*((v_new_ij dot e_ij) - (v_ij dot e_ij))
	     // i.e., only multiplication with vector r_ij is missing
	     newAbsRelVel *= 0.5 / pair -> absSquare();

// 	     MSG_DEBUG("ThermostatLA::thermalize", " REL VEL SETTING   rel vel = " << newAbsRelVel << endl << "pair->absSquare = " << pair->absSquare());
	     // now, delta[i] is computed and added to the old velocities
	     point_t temp;
	     temp = newAbsRelVel * pair->cartesian();
	     // it acts at least on one, right?
	     
// 	     if (pair->firstPart()->mySlot == 0)
// {
//   Particle* i = pair->firstPart();
//   MSG_DEBUG("ThermostatLA::thermalize", "BEFORE actsOnFirst" << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1]);
// }	   
// if (pair->secondPart()->mySlot == 0)
// {
//   Particle* i = pair->secondPart();
//   MSG_DEBUG("ThermostatLA::thermalize", "BEFORE actsOnFirst" << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1]);
// }
	     


	     
	     if (pair->actsOnFirst())
	       {
		 // here we have to preserve temp for the next if
		 point_t temp2;
		 m_firstAdd(&temp2, &(*pair));
		 for(size_t i = 0; i < SPACE_DIMS; ++i)
		   temp2[i]+=temp[i];
		 point_t temp3;
		 m_firstFac(&temp3, &(*pair));
		 for(size_t i = 0; i < SPACE_DIMS; ++i)
		   temp3[i]*=temp[i];
		 pair->firstPart()->v += temp3;
	       }
	     if (pair->actsOnSecond())
	       {
		 point_t temp2;
		 m_secondAdd(&temp2, &(*pair));
                 // now temp may be modified itself
		 for(size_t i = 0; i < SPACE_DIMS; ++i)
		   temp[i]*=temp2[i];
		 m_secondFac(&temp2, &(*pair));
		 for(size_t i = 0; i < SPACE_DIMS; ++i)
		   temp[i]*=temp2[i];
		 pair->secondPart()->v -= temp;
	       }
	       
	       
/*	       
	   if (pair->firstPart()->mySlot == 183)
{
  Particle* i = pair->firstPart();
  MSG_DEBUG("ThermostatLA::thermalize", "AFTER actsOnFirst" << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1]);
}	   
if (pair->secondPart()->mySlot == 183)
{
  Particle* i = pair->secondPart();
  MSG_DEBUG("ThermostatLA::thermalize", "AFTER actsOnFirst" << endl << "r=" << i->r << endl << "v=" << i->v << endl << "dt="  << i->dt << endl << "f0="  << i->force[0] << endl << "f1="  << i->force[1]);
}    */
	       
	       
	       
	       

#ifdef ENABLE_PTHREADS
	     pair->secondPart()->unlock();
	     pair->firstPart()->unlock();
#endif

//  MSG_DEBUG("ThermostatLA::thermalize", "AFTER THE LOOP!!!!!!!!!!  pair distance = " << pair->abs() << endl << "velocity  1st PART!!! = " << pair -> firstPart() -> v << endl << "slot  1st PART!!! = " << pair -> firstPart() -> mySlot << endl << "velocity  2nd PART!!! = " << pair -> secondPart() -> v << endl << "slot  2nd PART!!! = " << pair -> secondPart() -> mySlot);

	   }
       }
     );
}

