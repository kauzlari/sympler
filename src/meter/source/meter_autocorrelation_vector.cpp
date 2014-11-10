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


#include "meter_autocorrelation_vector.h"

#include "phase.h"
#include "simulation.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterAutocorrelationVector> meter_autocorrelation_vector("MeterAutocorrelationVector");

#define M_SIMULATION ((Simulation*)m_parent)
#define M_MANAGER M_SIMULATION->phase()->manager()


//---- Constructors/Destructor ----

MeterAutocorrelationVector::MeterAutocorrelationVector(Simulation *simulation)
  : Meter(simulation)
{                                               
  init();
}


MeterAutocorrelationVector::~MeterAutocorrelationVector()
{
}


void MeterAutocorrelationVector::init()
{
  m_properties.setClassName("MeterAutocorrelationVector");

  m_properties.setDescription
    ("This meter computes the autocorrelation funtion (ACF) for an arbitrary vector-property V stored in the Particle 's tag, i.e. <V.V(t)>, where the dot denotes the scalar product. Note that the autocorrelation of an average vector quantity of the WHOLE SYSTEM CANNOT be computed by this Meter.");

  STRINGPC
    (species,
     m_species,
     "Species containing the input data.");

  m_species = "UNDEF";

  STRINGPC
    (symbol,
     m_inputSymbol,
     "Symbol of the input value for which the ACF should be computed.");

  m_inputSymbol = "UNDEF";

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
  m_format.addAttribute("acf", DataFormat::DOUBLE);

}


//---- Methods ----

void MeterAutocorrelationVector::setup()
{
  Meter::setup();

  if(m_limitAcfAv <= 0)
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nOfACFs\"! Current value: " + ObjToString(m_limitAcfAv));

  if(m_nBuffAcf <= 0)
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nBuffs\"! Current value: " + ObjToString(m_nBuffAcf));

  if(m_nValAcf <= 0)
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nValAcf\"! Current value: " + ObjToString(m_nValAcf));

  if(m_nValAcf < m_nBuffAcf)
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Choose other value for attribute \"nBuffs\"! Current value: " + ObjToString(m_nBuffAcf) + ". Attribute \"nBuffs\" must be <= \"nValAcf\" (=" + ObjToString(m_nValAcf) + ").");

  int minSteps = m_from_step_on + int((m_measure_every_n*m_nValAcf)*(m_limitAcfAv+m_nBuffAcf-1)/m_nBuffAcf);
  if (M_SIMULATION->controller()->timesteps() < minSteps)
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Increase the number of timesteps (in the Controller) to at least " + ObjToString(minSteps) + "!");

  if (m_species == "UNDEF")
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "\"species\" undefined.");

  if (m_inputSymbol == "UNDEF")
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "\"symbol\" undefined.");

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

  if(!Particle::s_tag_format[m_colour].attrExists(m_inputSymbol))
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Symbol " + m_inputSymbol + " not found for species '" + M_MANAGER->species(m_colour) + "'!");
  if(Particle::s_tag_format[m_colour].attrByName(m_inputSymbol).datatype != DataFormat::POINT)
    throw gError("MeterAutocorrelationVector::setup"+FILE_INFO, "Symbol " + m_inputSymbol + " of species '" + M_MANAGER->species(m_colour) + "' is not a vector!");

  m_inputOffset = Particle::s_tag_format[m_colour].offsetByName(m_inputSymbol);

  string name = "MeterAutoCorrelation" ;
  name +=  "_" + m_inputSymbol + "_" + m_species;


  while(Particle::s_tag_format[m_colour].attrExists(name))
    name = "_" + name;

  m_acfOrgOffset = Particle::s_tag_format[m_colour].addAttribute(name, DataFormat::VECTOR_POINT, true, name).offset;

  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
    m_indexAcf[nb] = -int((nb)*m_nValAcf / m_nBuffAcf + 1);
//     MSG_DEBUG("MeterAutocorrelationVector::setup", "nb: " << nb << ": m_indexAcf[nb]=" << m_indexAcf[nb]);
  }
  zeroAcf();

}

void MeterAutocorrelationVector::zeroAcf()
{
  m_timeDone = false;
  m_countAcfAv = 0;
  for(size_t j = 0; j < m_nValAcf; ++j)
    m_acfAv[j] = 0.;
}


void MeterAutocorrelationVector::measureNow(const double& time)
{
  Phase *phase = ((Simulation*) m_parent)->phase();

  size_t nOfParticles = phase -> returnNofPartC(m_colour);
  
  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
//     MSG_DEBUG("MeterAutocorrelationVector::measureNow", "nb="<< nb << ", before incr: " << m_indexAcf[nb]);
    ++(m_indexAcf[nb]);
    if(m_indexAcf[nb] < 0) continue;
    if(m_indexAcf[nb] == 0) {
      FOR_EACH_PARTICLE_C 
	(phase, m_colour,

	 // FIXME: the following is not nice, but currently we can neither create the dynamic arrays immediately with the right size nor do we have non-dynamic arrays with fixed predefined size which are directly saved in the tag (not just the pointer)

	 __iSLFE->tag.vectorPointByOffset(m_acfOrgOffset)->resize(m_nBuffAcf);
	 
	 (*(__iSLFE->tag.vectorPointByOffset(m_acfOrgOffset)))[nb] = __iSLFE->tag.pointByOffset(m_inputOffset);

	 );
    }

    size_t ni = m_indexAcf[nb];
    if(nb == 0 && !m_timeDone) m_time[ni] = time;
    m_acf[nb][ni] = 0.;
//     MSG_DEBUG("MeterAutocorrelationVector::measureNow", "before second p-loop, nb=" << nb << ", ni=" << ni);
    FOR_EACH_PARTICLE_C 
      (phase, m_colour,
       // scalar product!
       m_acf[nb][ni] 
        += (*(__iSLFE->tag.vectorPointByOffset(m_acfOrgOffset)))[nb] * __iSLFE->tag.pointByOffset(m_inputOffset);

       );    

    // divide by number of particles
    m_acf[nb][ni] /= nOfParticles;
  }  
  accumAcf();
}


void MeterAutocorrelationVector::accumAcf(/*const double& time*/)
{
//   Phase *phase = ((Simulation*) m_parent)->phase();

  for(size_t nb = 0; nb < m_nBuffAcf; ++nb) {
    if(size_t(m_indexAcf[nb]) == m_nValAcf-1) {
      m_timeDone = true;
      for(size_t j = 0; j < m_nValAcf; ++j) {
	m_acfAv[j] += m_acf[nb][j];
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
	double firstVal = m_acfAv[0];
	double firstTime = m_time[0];
 	data->doubleByIndex(0) = 0.;
 	data->doubleByIndex(1) = firstVal / m_limitAcfAv/*1.*/;
 	distribute(data);
	for(size_t j = 1; j < m_nValAcf; ++j) {
 	  data->doubleByIndex(0) = m_time[j] - firstTime;
	  // divide by number of single CFs
	  data->doubleByIndex(1) = m_acfAv[j] / m_limitAcfAv; 
	  distribute(data);
	}

	zeroAcf();
      } // end: if(m_countAcfAv == m_limitAcfAv)
    } // end: if(size_t(m_indexAcf[nb]) == m_nValAcf-1)
  } // end: for(size_t nb = 0; nb < m_nBuffAcf; ++nb 
}
