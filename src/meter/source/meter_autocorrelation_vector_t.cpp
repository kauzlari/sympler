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


#include "meter_autocorrelation_vector_t.h"

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterAutocorrelationVectorT> meter_autocorrelation_vector_t("MeterAutocorrelationVectorT");

#define M_SIMULATION ((Simulation*)m_parent)
#define M_MANAGER M_SIMULATION->phase()->manager()


//---- Constructors/Destructor ----

MeterAutocorrelationVectorT::MeterAutocorrelationVectorT(Simulation *simulation)
  : Meter(simulation)
{                                               
  init();
}


MeterAutocorrelationVectorT::~MeterAutocorrelationVectorT()
{
}


void MeterAutocorrelationVectorT::init()
{
  m_properties.setClassName("MeterAutocorrelationVectorT");

  m_properties.setDescription
    ("This meter computes the tensorial autocorrelation funtion (ACF) for two arbitrary vector-properties V1(t=0) and V2(t) stored in the Particle's tag, i.e. the matrix\n"
"((<V1xV2x(t)>, <V1xV2y(t)>, <V1xV2z(t)>)\n"
" (<V1yV2x(t)>, <V1yV2y(t)>, <V1yV2z(t)>)\n"
" (<V1zV2x(t)>, <V1zV2y(t)>, <V1zV2z(t)>)).\n" 
"Note that the autocorrelation of an average vector quantity of the WHOLE SYSTEM CANNOT be computed by this Meter.");

  STRINGPC
    (species,
     m_species,
     "Species containing the input data.");

  m_species = "UNDEF";

  STRINGPC
    (symbol1,
     m_inputSymbol.first,
     "Symbol of the input value V1(t=0) for which the ACF should be computed.");

  m_inputSymbol.first = "UNDEF";

  STRINGPC
    (symbol2,
     m_inputSymbol.second,
     "Symbol of the input value V2(t) for which the ACF should be computed.");

  m_inputSymbol.second = "UNDEF";

  INTPC
    (nOfACFs, m_limitAcfAv, 0,
     "The number of individual ACFs to average over.");

  m_limitAcfAv = 0;

  INTPC
    (nBuffs, m_nBuffAcf, 0,
     "The number of buffers for simultaneously accumulated ACFs.");

  m_nBuffAcf = 0;

  INTPC
    (nValACF, m_nValAcf, 0,
     "The length of the ACF.");

  m_nValAcf = 0;



	// OLD
//   m_format.addAttribute("time", DataFormat::VECTOR_DOUBLE);
//   m_format.addAttribute("acf", DataFormat::VECTOR_DOUBLE);
  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute("acfxx", DataFormat::DOUBLE);
  m_format.addAttribute("acfxy", DataFormat::DOUBLE);
  m_format.addAttribute("acfxz", DataFormat::DOUBLE);
  m_format.addAttribute("acfyx", DataFormat::DOUBLE);
  m_format.addAttribute("acfyy", DataFormat::DOUBLE);
  m_format.addAttribute("acfyz", DataFormat::DOUBLE);
  m_format.addAttribute("acfzx", DataFormat::DOUBLE);
  m_format.addAttribute("acfzy", DataFormat::DOUBLE);
  m_format.addAttribute("acfzz", DataFormat::DOUBLE);

}


//---- Methods ----

void MeterAutocorrelationVectorT::setup()
{
  Meter::setup();

  if(m_limitAcfAv <= 0)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Choose other value for attribute \"nOfACFs\"! Current value: " + ObjToString(m_limitAcfAv));

  if(m_nBuffAcf <= 0)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Choose other value for attribute \"nBuffs\"! Current value: " + ObjToString(m_nBuffAcf));

  if(m_nValAcf <= 0)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Choose other value for attribute \"nValAcf\"! Current value: " + ObjToString(m_nValAcf));

  if(m_nValAcf < m_nBuffAcf)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Choose other value for attribute \"nBuffs\"! Current value: " + ObjToString(m_nBuffAcf) + ". Attribute \"nBuffs\" must be <= \"nValAcf\" (=" + ObjToString(m_nValAcf) + ").");

  int minSteps = m_from_step_on + int((m_measure_every_n*m_nValAcf)*(m_limitAcfAv+m_nBuffAcf-1)/m_nBuffAcf);
  if (M_SIMULATION->controller()->timesteps() < minSteps)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Increase the number of timesteps (in the Controller) to at least " + ObjToString(minSteps) + "!");

  if (m_species == "UNDEF")
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "\"species\" undefined.");

  if (m_inputSymbol.first == "UNDEF")
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "\"symbol1\" undefined.");

  if (m_inputSymbol.second == "UNDEF")
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "\"symbol2\" undefined.");

  m_colour = M_MANAGER->getColour(m_species);

  m_indexAcf.resize(m_nBuffAcf);

  m_acf.resize(m_nBuffAcf);
  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
    m_acf[nb].resize(m_nValAcf);
  }

  m_acfAv.resize(m_nValAcf);

  m_time.resize(m_nValAcf);

  // not needed because only one offset to a DataFormat::VECTOR_INT
//   m_acfOrgOffset.resize(m_nBuffAcf);

  if(!Particle::s_tag_format[m_colour].attrExists(m_inputSymbol.first))
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Symbol " + m_inputSymbol.first + " not found for species '" + M_MANAGER->species(m_colour) + "'!");
  if(Particle::s_tag_format[m_colour].attrByName(m_inputSymbol.first).datatype != DataFormat::POINT)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Symbol " + m_inputSymbol.first + " of species '" + M_MANAGER->species(m_colour) + "' is not a vector!");

  m_inputOffset.first = Particle::s_tag_format[m_colour].offsetByName(m_inputSymbol.first);

  if(!Particle::s_tag_format[m_colour].attrExists(m_inputSymbol.second))
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Symbol " + m_inputSymbol.second + " not found for species '" + M_MANAGER->species(m_colour) + "'!");
  if(Particle::s_tag_format[m_colour].attrByName(m_inputSymbol.second).datatype != DataFormat::POINT)
    throw gError("MeterAutocorrelationVectorT::setup"+FILE_INFO, "Symbol " + m_inputSymbol.second + " of species '" + M_MANAGER->species(m_colour) + "' is not a vector!");

  m_inputOffset.second = Particle::s_tag_format[m_colour].offsetByName(m_inputSymbol.second);

  string name = "MeterAutoCorrelation" ;
  name +=  "_" + m_inputSymbol.first + "_" + m_species;


  while(Particle::s_tag_format[m_colour].attrExists(name))
    name = "_" + name;

  m_acfOrgOffset = Particle::s_tag_format[m_colour].addAttribute(name, DataFormat::VECTOR_POINT, true, name).offset;

  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
    m_indexAcf[nb] = -int((nb)*m_nValAcf / m_nBuffAcf + 1);
//     MSG_DEBUG("MeterAutocorrelationVectorT::setup", "nb: " << nb << ": m_indexAcf[nb]=" << m_indexAcf[nb]);
  }
  zeroAcf();

}

void MeterAutocorrelationVectorT::zeroAcf()
{
  m_timeDone = false;
  m_countAcfAv = 0;

  for(size_t j = 0; j < m_nValAcf; ++j) {
    for(size_t k = 0; k < SPACE_DIMS; ++k) {
      for(size_t l = 0; l < SPACE_DIMS; ++l) {
	m_acfAv[j](k,l) = 0.;
      }
    }
  }

}


void MeterAutocorrelationVectorT::measureNow(const double& time)
{
  Phase *phase = ((Simulation*) m_parent)->phase();

  size_t nOfParticles = phase -> returnNofPartC(m_colour);
  
  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
//     MSG_DEBUG("MeterAutocorrelationVectorT::measureNow", "nb="<< nb << ", before incr: " << m_indexAcf[nb]);
    ++(m_indexAcf[nb]);
    if(m_indexAcf[nb] < 0) continue;
    if(m_indexAcf[nb] == 0) {
      FOR_EACH_PARTICLE_C 
	(phase, m_colour,

	 // FIXME: the following is not nice, but currently we can neither create the dynamic arrays immediately with the right size nor do we have non-dynamic arrays with fixed predefined size which are directly saved in the tag (not just the pointer)

	 __iSLFE->tag.vectorPointByOffset(m_acfOrgOffset)->resize(m_nBuffAcf);
	 
	 (*(__iSLFE->tag.vectorPointByOffset(m_acfOrgOffset)))[nb] = __iSLFE->tag.pointByOffset(m_inputOffset.first);

	 );
    }

    size_t ni = m_indexAcf[nb];
    if(nb == 0 && !m_timeDone) m_time[ni] = time;

    for(size_t k = 0; k < SPACE_DIMS; ++k) {
      for(size_t l = 0; l < SPACE_DIMS; ++l) {
	m_acf[nb][ni](k,l) = 0.;
      }
    }

    FOR_EACH_PARTICLE_C 
      (phase, m_colour,

       for(size_t k = 0; k < SPACE_DIMS; ++k) {
	 for(size_t l = 0; l < SPACE_DIMS; ++l) {
	   m_acf[nb][ni](k,l)
	     +=
	     // first term: kth component of origin of particle i
	     ((*(__iSLFE->tag.vectorPointByOffset(m_acfOrgOffset)))[nb])[k] 
	     // second term: lth component of current value of particle j
	     * (__iSLFE->tag.pointByOffset(m_inputOffset.second))[l];
	   
	 }
       }   

       );    

    // divide by number of particles
    m_acf[nb][ni] /= nOfParticles;
  }  
  accumAcf();
}


void MeterAutocorrelationVectorT::accumAcf(/*const double& time*/)
{
//   Phase *phase = ((Simulation*) m_parent)->phase();

  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
    if(size_t(m_indexAcf[nb]) == m_nValAcf-1) {
      m_timeDone = true;

      for(size_t j = 0; j < m_nValAcf; ++j) {
	for(size_t k = 0; k < SPACE_DIMS; ++k) {
	  for(size_t l = 0; l < SPACE_DIMS; ++l) {
	    (m_acfAv[j])(k,l) += (m_acf[nb][j])(k,l);
	  }
	}
      }

      m_indexAcf[nb] = -1;
      ++m_countAcfAv;
      if(m_countAcfAv == m_limitAcfAv) {

	// the integral (TAKE CARE: the result must still be multiplied by the timestep dt !!! Currently it is assumed that this is done EXTERNALY !)
// 	double fac = 1. / (SPACE_DIMS*phase->nParticles(m_colour)*m_limitAcfAv);
// 	m_acfInt = fac * m_stepAcf * integrate(m_acfAv, m_nValAcf);

	// store the NON-normalised acf
	// NEW: we distribute the data as single doubles so that we can get nice x and y columns
	data_sp data = m_format.newData();
	double firstTime = m_time[0];

	for(size_t j = 0; j < m_nValAcf; ++j) {
 	  data->doubleByIndex(0) = m_time[j] - firstTime;
	  // divide by number of single CFs
	  data->doubleByIndex(1) = m_acfAv[j](0,0) / m_limitAcfAv; 
	  data->doubleByIndex(2) = m_acfAv[j](0,1) / m_limitAcfAv; 
	  data->doubleByIndex(3) = m_acfAv[j](0,2) / m_limitAcfAv; 
	  data->doubleByIndex(4) = m_acfAv[j](1,0) / m_limitAcfAv; 
	  data->doubleByIndex(5) = m_acfAv[j](1,1) / m_limitAcfAv; 
	  data->doubleByIndex(6) = m_acfAv[j](1,2) / m_limitAcfAv; 
	  data->doubleByIndex(7) = m_acfAv[j](2,0) / m_limitAcfAv; 
	  data->doubleByIndex(8) = m_acfAv[j](2,1) / m_limitAcfAv; 
	  data->doubleByIndex(9) = m_acfAv[j](2,2) / m_limitAcfAv; 
	  distribute(data);
	}

	zeroAcf();
      } // end: if(m_countAcfAv == m_limitAcfAv)
    } // end: if(size_t(m_indexAcf[nb]) == m_nValAcf-1)
  } // end: for(size_t nb = 0; nb < m_nBuffAcf; ++nb 
}
