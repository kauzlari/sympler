/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


extern "C" {
  #include <freesteam/steam_pT.h>
}
#include "density_calculation.h"
#include "simulation.h"
#include "manager_cell.h"
#include <iostream>
#include <string>

const SymbolRegister<DensityCalculation> DensityCalculation("DensityCalculation");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

DensityCalculation::DensityCalculation
  (size_t colour, size_t offset, string symbolName)
  : ParticleCache(colour, offset, symbolName)/*m_colour(colour), m_offset(offset), *//*m_wf(wf)*/ {
  m_datatype = DataFormat::DOUBLE;
}

DensityCalculation::DensityCalculation
    (/*Node*/Simulation* parent)
  : ParticleCache(parent) {
  m_datatype = DataFormat::DOUBLE;
  init();
}


DensityCalculation::~DensityCalculation() {
  if(m_array_rho) {
    for (int i = 0; i < m_arraysize_pressure; ++i) {
      if(m_array_rho[i])
	delete [] m_array_rho[i];
    }
    delete [] m_array_rho;
  }
}

void DensityCalculation::setupLUT(double Tmin, double pmin, double Tmax, double pmax, int m_arraysize_pressure, int m_arraysize_temperature) {
  // Function intern auxiliary vaiables
  double pr = 0;
  double t = 0;
  double density;
  // Step sizes depending on the size of the Array and the ranges of the input values.
  m_calcstepP= (pmax-pmin)/(m_arraysize_pressure-1);
  m_calcstepT= (Tmax-Tmin)/(m_arraysize_temperature-1);
  // Initialization of the 2D Array, depending on the given temperature and pressure sizes.
  m_array_rho = new double*[m_arraysize_pressure+1];
  for (int i = 0; i <= m_arraysize_pressure; i++) {
    m_array_rho[i] = new double [m_arraysize_temperature+1];
  }
  // Loop calculates the densities in an adjusted temperature and pressure range and
  // stores them into the LUT.
  for(int j = 0; j < m_arraysize_pressure; j++) {
    for (int i = 0; i < m_arraysize_temperature; i++) {  
      SteamState S = freesteam_set_pT(pmin + pr,Tmin + t);
      density = freesteam_rho(S);
      m_array_rho[j][i] = density;
      t += m_calcstepT;
    }
    t = 0;
    pr += m_calcstepP;
  }
}

double DensityCalculation::calculateDensity(double inputT, double inputP, double Tmin, double pmin) {
  // Calculation of the surrounded sampling points    .
  int x_pressure_array_0 = (floor((inputP-pmin)/m_calcstepP));
  int y_temperature_array_0 =(floor((inputT-Tmin)/m_calcstepT));
  int x_pressure_array_1 = (floor((inputP-pmin)/m_calcstepP)+1);
  int y_temperature_array_1 =(floor((inputT-Tmin)/m_calcstepT)+1);
  // Normalization of the input values.
  double press_nominated =  ((inputP -(x_pressure_array_0*m_calcstepP+pmin)) / m_calcstepP);
  double temp_nominated =   ((inputT -(y_temperature_array_0*m_calcstepT+Tmin)) / m_calcstepT);
  // Bilinear interpolation of the normalized values.
  double x_0_y_0 = m_array_rho[x_pressure_array_0][y_temperature_array_0]*(1-press_nominated)*(1-temp_nominated);
  double x_1_y_0 = m_array_rho[x_pressure_array_1][y_temperature_array_0]*press_nominated*(1-temp_nominated);
  double x_0_y_1 = m_array_rho[x_pressure_array_0][y_temperature_array_1]*(1-press_nominated)*temp_nominated;
  double x_1_y_1 = m_array_rho[x_pressure_array_1][y_temperature_array_1]*temp_nominated*press_nominated;
  // Interpolated density value.
  m_density_interpolation = x_0_y_0 + x_1_y_0 + x_0_y_1 + x_1_y_1;
  return m_density_interpolation;
}



void DensityCalculation::init() {
  m_properties.setClassName("DensityCalculation");
  m_properties.setName("DensityCalculation");
  m_properties.setDescription("Contribution of the particle itself to its local density.");  
  STRINGPC
      (symbol, m_symbolName,
       "Name for the local density.");
  STRINGPC
      (pressure, m_pressureName,
       "Name for the local pressure.");
  STRINGPC
      (temperature, m_temperatureName,
       "Name for the local temperature.");
  DOUBLEPC
      (temperatureMax, m_Tmax, 0,
       "The local maximum temperature.");
  DOUBLEPC
      (temperatureMin, m_Tmin, 0,
       "The local maximum temperature.");
  DOUBLEPC
      (pressureMax, m_pmax, 0,
       "The local maximum pressure.");
  DOUBLEPC
      (pressureMin, m_pmin, 0,
       "The local minimum pressure.");
  INTPC
      (arraysize_pressure, m_arraysize_pressure, 0,
       "The size of the array for pressure values.");
  INTPC
      (arraysize_temperature, m_arraysize_temperature, 0,
       "The size of the array for temperature values.");

  m_temperatureName = "undefined";
  m_pressureName = "undefined";
  m_Tmin = 0;
  m_pmin = 0;
  m_Tmax = 0;
  m_pmax = 0;
  m_arraysize_pressure = 0;
  m_arraysize_temperature = 0;
}

void DensityCalculation::setup() {
  
  // Checks if all necessary input values are defined.
  if(m_species == "undefined")
    throw gError("DensityCalculation::setup", "Attribute 'species' has value \"undefined\""); 
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("DensityCalculation::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_temperatureName == "undefined")
    throw gError("DensityCalculation::setup", "Attribute 'temperature' has value \"undefined\""); 
  if(m_pressureName == "undefined")
    throw gError("DensityCalculation::setup", "Attribute 'pressure' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_Tmin == 0)
    throw gError("DensityCalculation::setup", "Attribute 'temperatureMin' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_Tmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'temperatureMax' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmin == 0)
    throw gError("DensityCalculation::setup", "Attribute 'pressureMin' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'pressureMax' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'arraysize_temperature' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_pmax == 0)
    throw gError("DensityCalculation::setup", "Attribute 'arraysize_pressure' has none of the allowed values \"0\", \"1\", \"2\".");

  
  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours, too?
  if(m_species == "ALL")
  {
    // YES! Using m_colour in this loop is correct since we want to
    // create one ParticleCache per colour
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      {
	throw gError("DensityCalculation::setup", "Symbol '" + m_symbolName + "' was already created by other module.");
      }
      else
        m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;

      if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName)) {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
          throw gError("DensityCalculation::setup", "Symbol '" + m_temperatureName + "' already exists as a non-scalar.");
        else m_temperatureOffset = Particle::s_tag_format[m_colour].offsetByName(m_temperatureName);
      }
      else 
	throw gError("DensityCalculation::setup", "Symbol '" + m_temperatureName + "' does not exist but required by this module.");
     
      if(Particle::s_tag_format[m_colour].attrExists(m_pressureName)) {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_pressureName).datatype)
          throw gError("DensityCalculation::setup", "Symbol '" + m_pressureName + "' already exists as a non-scalar.");
        else m_pressureOffset = Particle::s_tag_format[m_colour].offsetByName(m_pressureName);
      }
      else 
	throw gError("DensityCalculation::setup", "Symbol '" + m_pressureName + "' does not exist but required by this module.");

      // is it the last cache to be created?
      if(m_colour == M_MANAGER->nColours()-1)
      {
        if(m_phaseUser == 0)
          Particle::registerCache_0(this);
        else if(m_phaseUser == 1) 
          Particle::registerCache(this);
        else // so it is 2
        {
          ParticleCache* pc = new DensityCalculation(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((DensityCalculation*) pc)->stage() == m_stage);
          assert(((DensityCalculation*) pc)->m_colour == m_colour);
          assert(((DensityCalculation*) pc)->m_offset == m_offset);
	  assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	  assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);
          Particle::registerCache(pc);
          Particle::registerCache_0(this);
        }
      }
      // No? Then make a copy
      else 
      {
        ParticleCache* pc = new DensityCalculation(*this);
        assert(pc->mySymbolName() == m_symbolName);
        assert(((DensityCalculation*) pc)->stage() == m_stage);
        assert(((DensityCalculation*) pc)->m_colour == m_colour);
        assert(((DensityCalculation*) pc)->m_offset == m_offset);
	assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);

        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);
          
          pc = new DensityCalculation(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((DensityCalculation*) pc)->stage() == m_stage);
          assert(((DensityCalculation*) pc)->m_colour == m_colour);
          assert(((DensityCalculation*) pc)->m_offset == m_offset);
	  assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	  assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);
          
          Particle::registerCache_0(pc);
        }
      
      }
    }
  }
  else /*it is a Symbol limited to one colour*/
  {
    m_colour = M_MANAGER->getColour(m_species);
    
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
    {
      throw gError("ParticleCacheDensitySelfContribution::setup", "Symbol '" + m_symbolName + "' was already created by other module.");
    }
    else
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;

    if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName)) {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
        throw gError("PressureCalculation::setup", "Symbol '" + m_temperatureName + "' already exists as a non-scalar.");
      else m_temperatureOffset = Particle::s_tag_format[m_colour].offsetByName(m_temperatureName);
    } else 
      throw gError("PressureCalculation::setup", "Symbol '" + m_temperatureName + "' does not exist but required by this module.");

    if(Particle::s_tag_format[m_colour].attrExists(m_pressureName)) {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_pressureName).datatype)
	throw gError("PressureCalculation::setup", "Symbol '" + m_pressureName + "' already exists as a non-scalar.");
      else m_pressureOffset = Particle::s_tag_format[m_colour].offsetByName(m_pressureName);
    } else 
      throw gError("PressureCalculation::setup", "Symbol '" + m_pressureName + "' does not exist but required by this module.");

    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else
    {
      ParticleCache* pc = new DensityCalculation(*this);
      assert(pc->mySymbolName() == m_symbolName);
      assert(((DensityCalculation*) pc)->stage() == m_stage);
      assert(((DensityCalculation*) pc)->m_colour == m_colour);
      assert(((DensityCalculation*) pc)->m_offset == m_offset);
      assert(((DensityCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
      assert(((DensityCalculation*) pc)->m_pressureOffset == m_pressureOffset);

      Particle::registerCache(pc);
      Particle::registerCache_0(this);
    } 
  }

   setupLUT( m_Tmin,  m_pmin, m_Tmax, m_pmax, m_arraysize_pressure, m_arraysize_temperature);   
}


void DensityCalculation::registerWithParticle()
{
}

