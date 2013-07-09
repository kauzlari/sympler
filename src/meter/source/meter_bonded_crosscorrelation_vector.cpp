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


#include "meter_bonded_crosscorrelation_vector.h"

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterBondedCrosscorrelationVector> meter_bonded_crosscorrelation_vector("MeterBondedCrosscorrelationVector");

#define M_SIMULATION ((Simulation*)m_parent)
#define M_MANAGER M_SIMULATION->phase()->manager()


//---- Constructors/Destructor ----

MeterBondedCrosscorrelationVector::MeterBondedCrosscorrelationVector(Simulation *simulation)
  : Meter(simulation)
{                                               
  init();
}


MeterBondedCrosscorrelationVector::~MeterBondedCrosscorrelationVector()
{
}


void MeterBondedCrosscorrelationVector::init()
{
  m_properties.setClassName("MeterBondedCrosscorrelationVector");

  m_properties.setDescription
    ("This meter computes the crosscorrelation funtion (CCF) for two arbitrary vector-properties V1(t=0), V2(t) defined by the attributes \"symbol1\" and \"symbol2\" and stored in the tag of a pair of Particles P1, P2 whose species are defined by the attributes \"species1\" and \"species2\". The outer product of the vectors is computed resulting in a tensor. The two particles belong to a list of bonded pairs specified by the user.");

  STRINGPC(listName, m_listName, "Identifier of the list of bonded pairs, this Meter should work on.");

  m_listName = "undefined";

  STRINGPC
    (species1,
     m_species.first,
     "Species containing the input data.");

  STRINGPC
    (species2,
     m_species.second,
     "Species containing the input data.");

  m_species.first = "UNDEF";
  m_species.second = "UNDEF";

  STRINGPC
    (symbol1,
     m_inputSymbol.first,
     "Symbol of the input value V1(t=0) of particles of species1 for which the CCF should be computed.");

  STRINGPC
    (symbol2,
     m_inputSymbol.second,
     "Symbol of the input value V2(t) of particles of species2 for which the CCF should be computed.");

  m_inputSymbol.first = "UNDEF";
  m_inputSymbol.second = "UNDEF";

  INTPC
    (nOfCFs, m_limitCfAv, 0,
     "The number of individual CCFs to average over.");

  m_limitCfAv = 0;

  INTPC
    (nBuffs, m_nBuffCf, 0,
     "The number of buffers for simultaneously accumulated CCFs.");

  m_nBuffCf = 0;

  INTPC
    (nValCF, m_nValCf, 0,
     "The length of the CCF.");

  m_nValCf = 0;

  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute("cfxx", DataFormat::DOUBLE);
  m_format.addAttribute("cfxy", DataFormat::DOUBLE);
  m_format.addAttribute("cfxz", DataFormat::DOUBLE);
  m_format.addAttribute("cfyx", DataFormat::DOUBLE);
  m_format.addAttribute("cfyy", DataFormat::DOUBLE);
  m_format.addAttribute("cfyz", DataFormat::DOUBLE);
  m_format.addAttribute("cfzx", DataFormat::DOUBLE);
  m_format.addAttribute("cfzy", DataFormat::DOUBLE);
  m_format.addAttribute("cfzz", DataFormat::DOUBLE);

}


//---- Methods ----

void MeterBondedCrosscorrelationVector::setup()
{
  Meter::setup();


  if(m_listName == "undefined")
    throw gError("MeterBondedCrosscorrelationVector::setup", "Attribute 'listName' is undefined!");


  m_colours.first = M_MANAGER->getColour(m_species.first);
  m_colours.second = M_MANAGER->getColour(m_species.second);

  m_cp = M_MANAGER->cp(m_colours.first, m_colours.second);

  // should work here because the lists have been declared in setup() of modules creating connections (e.g., in pc_connector.cpp)
  m_listIndex = m_cp->connectedListIndex(m_listName);

  if(m_species.first == "undefined")
    throw gError("MeterBondedCrosscorrelationVector::setup", ": Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled."); 
  if(m_species.second == "undefined")
    throw gError("MeterBondedCrosscorrelationVector::setup", ": Attribute 'species1' has value \"undefined\" and 'allPairs' is disabled."); 



  if(m_limitCfAv <= 0)
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nOfCFs\"! Current value: " + ObjToString(m_limitCfAv));

  if(m_nBuffCf <= 0)
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nBuffs\"! Current value: " + ObjToString(m_nBuffCf));

  if(m_nValCf <= 0)
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nValCf\"! Current value: " + ObjToString(m_nValCf));

  if(m_nValCf < m_nBuffCf)
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nBuffs\"! Current value: " + ObjToString(m_nBuffCf) + ". Attribute \"nBuffs\" must be <= \"nValCf\" (=" + ObjToString(m_nValCf) + ").");

  int minSteps = m_from_step_on + int((m_measure_every_n*m_nValCf)*(m_limitCfAv+m_nBuffCf-1)/m_nBuffCf);
  if (M_SIMULATION->controller()->timesteps() < minSteps)
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Increase the number of timesteps (in the Controller) to at least " + ObjToString(minSteps) + "!");

  if (m_species.first == "UNDEF")
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "\"species1\" undefined.");

  if (m_species.second == "UNDEF")
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "\"species2\" undefined.");

  if (m_inputSymbol.first == "UNDEF")
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "\"symbol1\" undefined.");

  if (m_inputSymbol.second == "UNDEF")
    throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "\"symbol2\" undefined.");


  m_indexCf.resize(m_nBuffCf);

  m_cf.resize(m_nBuffCf);
  for(size_t nb = 0; nb < m_nBuffCf; ++nb) {
    m_cf[nb].resize(m_nValCf);
  }

  m_cfAv.resize(m_nValCf);

  m_time.resize(m_nValCf);

  // not needed because only one offset to a DataFormat::VECTOR_INT
//   m_cfOrgOffset.resize(m_nBuffCf);





 if(!Particle::s_tag_format[m_cp->firstColour()].attrExists(m_inputSymbol.first))
   throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Symbol " + m_inputSymbol.first + " not found for species '" + m_species.first + "'!");
 if(Particle::s_tag_format[m_cp->firstColour()].attrByName(m_inputSymbol.first).datatype != DataFormat::POINT)
   throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Symbol " + m_inputSymbol.first + " of species '" + m_species.first + "' is not a vector!");
 m_inputOffset.first = Particle::s_tag_format[m_cp->firstColour()].offsetByName(m_inputSymbol.first);

 if(!Particle::s_tag_format[m_cp->secondColour()].attrExists(m_inputSymbol.second))
   throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Symbol " + m_inputSymbol.second + " not found for species '" + m_species.second + "'!");
 if(Particle::s_tag_format[m_cp->secondColour()].attrByName(m_inputSymbol.second).datatype != DataFormat::POINT)
   throw gError("MeterBondedCrosscorrelationVector::setup"+FILE_INFO, "Symbol " + m_inputSymbol.second + " of species '" + m_species.second + "' is not a vector!");
 m_inputOffset.second = Particle::s_tag_format[m_cp->secondColour()].offsetByName(m_inputSymbol.second);


  string name = "MeterBondedCrossCorrelation";
  name +=  "_" + m_inputSymbol.first + "_" + m_species.first;

  while(Particle::s_tag_format[m_cp->firstColour()].attrExists(name))
    name = "_" + name;

  m_cfOrgOffset/*.first*/ = Particle::s_tag_format[m_cp->firstColour()].addAttribute(name, DataFormat::VECTOR_POINT, true, name).offset;

//---   next currently not necessary ---

//   // this must be repeated to get analogous names
//   name = "MeterBondedCrossCorrelation";
//   name +=  "_" + m_inputSymbol + "_" + m_species.second;

//   while(Particle::s_tag_format[m_cp->secondColour()].attrExists(name))
//     name = "_" + name;

//   m_cfOrgOffset.second = Particle::s_tag_format[m_cp->secondColour()].addAttribute(name, DataFormat::VECTOR_POINT, true, name).offset;

//END: ---   next currently not necessary ---


  for(size_t nb = 0; nb < m_nBuffCf; ++nb) {
    m_indexCf[nb] = -int((nb)*m_nValCf / m_nBuffCf + 1);
//     MSG_DEBUG("MeterBondedCrosscorrelationVector::setup", "nb: " << nb << ": m_indexCf[nb]=" << m_indexCf[nb]);
  }
  zeroCf();

}

void MeterBondedCrosscorrelationVector::zeroCf()
{
  m_timeDone = false;
  m_countCfAv = 0;
  for(size_t j = 0; j < m_nValCf; ++j) {
    for(size_t k = 0; k < SPACE_DIMS; ++k) {
      for(size_t l = 0; l < SPACE_DIMS; ++l) {
	m_cfAv[j](k,l) = 0.;
      }
    }
  }
}

void MeterBondedCrosscorrelationVector::measureNow(const double& time)
{
  Phase *phase = ((Simulation*) m_parent)->phase();

  PairList* myPairList = m_cp->connectedList(m_listIndex);

  size_t nOfPairs = myPairList -> size();
  
  for(size_t nb = 0; nb < m_nBuffCf; ++nb) {
    ++(m_indexCf[nb]);
    if(m_indexCf[nb] < 0) continue;
    if(m_indexCf[nb] == 0) {

      FOR_EACH_PARTICLE_C 
	(phase, m_colours.first,

	 __iSLFE->tag.vectorPointByOffset(m_cfOrgOffset)->resize(m_nBuffCf);
	 
	 (*(__iSLFE->tag.vectorPointByOffset(m_cfOrgOffset))/*tmpVec*/)[nb] = __iSLFE->tag.pointByOffset(m_inputOffset.first);
	 );
    }

    size_t ni = m_indexCf[nb];
    if(nb == 0 && !m_timeDone) m_time[ni] = time;
    for(size_t k = 0; k < SPACE_DIMS; ++k) {
      for(size_t l = 0; l < SPACE_DIMS; ++l) {
	m_cf[nb][ni](k,l) = 0.;
      }
    }
//     MSG_DEBUG("MeterBondedCrosscorrelationVector::measureNow", "before second p-loop, nb=" << nb << ", ni=" << ni);
    
    // loop over pairs of current connected list 
    for(Pairdist *pair = myPairList->first(); pair != NULL; pair = pair->next) {
      Particle* i = pair->firstPart();   
      Particle* j = pair->secondPart();   
      // dyadic product
      for(size_t k = 0; k < SPACE_DIMS; ++k) {
	for(size_t l = 0; l < SPACE_DIMS; ++l) {
	  m_cf[nb][ni](k,l)
	    +=
	    // first term: kth component of origin of particle i
	    ((*(i->tag.vectorPointByOffset(m_cfOrgOffset/*.first*/)))[nb])[k] 
	    // second term: lth component of current value of particle j
 	    * (j->tag.pointByOffset(m_inputOffset.second))[l];

	}
      }   
    }

    // divide by number of pairs
    m_cf[nb][ni] /= nOfPairs;
    
  } // end of: for(size_t nb = 0; nb < m_nBuffCf; ++nb) 
  accumCf();
}


void MeterBondedCrosscorrelationVector::accumCf(/*const double& time*/)
{
//   Phase *phase = ((Simulation*) m_parent)->phase();

  for(size_t nb = 0; nb < m_nBuffCf; ++nb) {
    if(size_t(m_indexCf[nb]) == m_nValCf-1) {
      m_timeDone = true;

      for(size_t j = 0; j < m_nValCf; ++j) {
	for(size_t k = 0; k < SPACE_DIMS; ++k) {
	  for(size_t l = 0; l < SPACE_DIMS; ++l) {
	    (m_cfAv[j])(k,l) += (m_cf[nb][j])(k,l);
	  }
	}
      }

      m_indexCf[nb] = -1;
      ++m_countCfAv;
      if(m_countCfAv == m_limitCfAv) {

	// store the NON-normalised cf
	// NEW: we distribute the data as single doubles so that we can get nice x and y columns
	data_sp data = m_format.newData();
	double firstTime = m_time[0];
	for(size_t j = 0; j < m_nValCf; ++j) {
 	  data->doubleByIndex(0) = m_time[j] - firstTime;
	  // divide by number of single CFs
	  data->doubleByIndex(1) = m_cfAv[j](0,0) / m_limitCfAv; 
	  data->doubleByIndex(2) = m_cfAv[j](0,1) / m_limitCfAv; 
	  data->doubleByIndex(3) = m_cfAv[j](0,2) / m_limitCfAv; 
	  data->doubleByIndex(4) = m_cfAv[j](1,0) / m_limitCfAv; 
	  data->doubleByIndex(5) = m_cfAv[j](1,1) / m_limitCfAv; 
	  data->doubleByIndex(6) = m_cfAv[j](1,2) / m_limitCfAv; 
	  data->doubleByIndex(7) = m_cfAv[j](2,0) / m_limitCfAv; 
	  data->doubleByIndex(8) = m_cfAv[j](2,1) / m_limitCfAv; 
	  data->doubleByIndex(9) = m_cfAv[j](2,2) / m_limitCfAv; 
	  distribute(data);
	}

	zeroCf();
      }
    }
  }
}


// next should not be necessary and work fine in setup()
// MeterBondedCrosscorrelationVector::setupAfterParticleCreation()
// {
//  
//   m_listIndex = m_cp->connectedListIndex(m_listName);

// }
