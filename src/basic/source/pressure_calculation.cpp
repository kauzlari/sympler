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


extern "C" {
  #include <freesteam/b23.h>
  #include <freesteam/steam_pT.h>
  #include <freesteam/region3.h>
}

#include "pressure_calculation.h"
#include "simulation.h"
#include "manager_cell.h"

const SymbolRegister<PressureCalculation> PressureCalculation("PressureCalculation");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PressureCalculation::PressureCalculation
  (size_t colour, size_t offset, string symbolName)
  : ParticleCache(colour, offset, symbolName) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
}

PressureCalculation::PressureCalculation
    (/*Node*/Simulation* parent)
  : ParticleCache(parent) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  init();
}


PressureCalculation::~PressureCalculation()
{
  if(m_array_p) {
    for (int i = 0; i < m_arraysize_density; ++i) {
      if(m_array_p[i])
	delete [] m_array_p[i];
    }
    delete [] m_array_p;
  }
}

void PressureCalculation::setupLUT() {

  // auxiliary vaiables
  long double rho = 0;
  long double deltaT = 0;
  double pressure;

  // Temperature boarders for Region 3. See IAPWS-IF97 for more information.
  if (m_Tmin > 623.15 && m_Tmax < 863.15) {
    // Step sizes depending on the size of the Array and the ranges of the input values.
    m_calcstepRho = (m_rhomax - m_rhomin)/(m_arraysize_density-1);
    m_calcstepT = (m_Tmax -m_Tmin)/(m_arraysize_temperature-1);
    // Initialization of the 2D Array, depending on the given temperature and pressure sizes.
    m_array_p = new double*[m_arraysize_density+1];
    for (int i = 0; i <= m_arraysize_density; i++) {
      m_array_p[i] = new double [m_arraysize_temperature+1];
    }
    // Minimum density boarder along b23-line.
    double densityBoundary;
    // Loop calculates the pressures in an adjusted temperature and density range and
    // stores them into the LUT.
    for(int j = 0; j < m_arraysize_density; j++) {
      for (int i = 0; i < m_arraysize_temperature; i++) {
        // Calculates the minimum and maximum density boarders to check
        // if the given density is in Region 3.	  
        double b23Pressure = freesteam_b23_p_T(deltaT + m_Tmin);
        SteamState S = freesteam_set_pT(b23Pressure, deltaT + m_Tmin);
        densityBoundary = freesteam_rho(S);
        if (m_rhomin + rho > densityBoundary) {
	  pressure = freesteam_region3_p_rhoT(m_rhomin + rho, m_Tmin + deltaT);
          m_array_p[j][i] = pressure;
          deltaT += m_calcstepT;
        } else {
	  throw gError("PressureCalculation::setup", "Requested density and temperature parameters aren't within region 3. (See IAPWS-IF97 for more information): Unsuccessfully tried to compute P(T=" + ObjToString(deltaT+m_Tmin) + ",rho=" + ObjToString(m_rhomin+rho) + ") crossing rho_min(T) = " + ObjToString(densityBoundary));
        }
      }
      deltaT=0;
      rho += m_calcstepRho;
    }
 
  } else {
     throw gError("PressureCalculation::setup", "Requested temperature parameters aren't within region 3 (623.15K < T < 863.15K). (See IAPWS-IF97 for more information)");
  }     
}

double PressureCalculation::calculatePressure(double inputT, double inputRho) {

  // out of bounds?
  if((inputRho < m_rhomin) || (inputT < m_Tmin) || (inputRho > m_rhomax) || (inputT > m_Tmax))
    throw gError("DensityCalculation::calculatePressure", "(rho,T) pair out of bounds: rho=" + ObjToString(inputRho) + ", T=" + ObjToString(inputT) + ", admissible [rhomin,rhomax] = [" + ObjToString(m_rhomin) + "," + ObjToString(m_rhomax) + "], admissible [Tmin,Tmax] = [" + ObjToString(m_Tmin) + "," + ObjToString(m_Tmax) + "].");

  // Calculation of the surrounding sampling points 
  int x_pressure_array_0 = (floor((inputRho-m_rhomin)/m_calcstepRho));
  int y_temperature_array_0 =(floor((inputT-m_Tmin)/m_calcstepT));
  int x_pressure_array_1 = (floor((inputRho-m_rhomin)/m_calcstepRho)+1);
  int y_temperature_array_1 =(floor((inputT-m_Tmin)/m_calcstepT)+1);
  // Normalisation of the input values.
  double dens_normalised =  ((inputRho - (x_pressure_array_0*m_calcstepRho + m_rhomin)) /m_calcstepRho);
  double temp_normalised =   ((inputT - (y_temperature_array_0*m_calcstepT + m_Tmin)) / m_calcstepT);
  // Bilinear interpolation of the normalized values.
  double x_0_y_0 = m_array_p[x_pressure_array_0][y_temperature_array_0]*(1-dens_normalised)*(1-temp_normalised);
  double x_1_y_0 = m_array_p[x_pressure_array_1][y_temperature_array_0]*dens_normalised*(1-temp_normalised);
  double x_0_y_1 = m_array_p[x_pressure_array_0][y_temperature_array_1]*(1-dens_normalised)*temp_normalised;
  double x_1_y_1 = m_array_p[x_pressure_array_1][y_temperature_array_1]*temp_normalised*dens_normalised;

  // Interpolated density value.
  return x_0_y_0 + x_1_y_0 + x_0_y_1 + x_1_y_1;

}

void PressureCalculation::init()
{
  m_properties.setClassName("PressureCalculation");
  m_properties.setName("PressureCalculation");
  m_properties.setDescription
    ("Local pressure at the particle computed from the local density "
     "and local temperature based on IAPWS-IF97 (International "
     "Association for the Properties of Water and Steam. Revised "
     "release on the IAPWS industrial formulation 1997 for the "
     "thermodynamic properties of water and steam. adadad, August "
     "2007)."
     ); 

  STRINGPC
      (symbol, m_symbolName,
       "Name for the computed local pressure.");
  STRINGPC
      (density, m_densityName,
       "Name for the local density used as input.");
  STRINGPC
      (temperature, m_temperatureName,
       "Name for the local temperature used as input.");
  DOUBLEPC
      (temperatureMax, m_Tmax, 0,
       "Upper limit for the admissible temperature range.");
  DOUBLEPC
      (temperatureMin, m_Tmin, 0,
       "Lower limit for the admissible temperature range.");
  DOUBLEPC
      (densityMax, m_rhomax, 0,
       "Upper limit for the admissible density range.");
  DOUBLEPC
      (densityMin, m_rhomin, 0,
       "Lower limit for the admissible density range.");
  INTPC
      (arraysize_density, m_arraysize_density, 0,
       "The size of the array for pressure values in the generated "
       "look-up table.");
  INTPC
      (arraysize_temperature, m_arraysize_temperature, 0,
       "The size of the array for temperature values in the generated "
       "look-up table.");
  
  m_temperatureName = "undefined";
  m_densityName = "undefined";
  m_Tmin = HUGE_VAL;
  m_rhomin = HUGE_VAL;
  m_Tmax = HUGE_VAL;
  m_rhomax = HUGE_VAL;
  m_arraysize_density = 0;
  m_arraysize_temperature = 0;
}

void PressureCalculation::setup()
{
  if(m_species == "undefined")
    throw gError("PressureCalculation::setup", "Attribute 'species' has value \"undefined\""); 
  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("PressureCalculation::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  if(m_temperatureName == "undefined")
    throw gError("PressureCalculation::setup", "Attribute 'temperature' has value \"undefined\""); 
  if(m_densityName == "undefined")
    throw gError("PressureCalculation::setup", "Attribute 'density' has value \"undefined\"");
  if(m_Tmin == HUGE_VAL)
    throw gError("PressureCalculation::setup", "Attribute 'temperatureMin' was not defined");
  if(m_Tmax == HUGE_VAL)
    throw gError("PressureCalculation::setup", "Attribute 'temperatureMax' was not defined");
  if(m_rhomin == HUGE_VAL)
    throw gError("PressureCalculation::setup", "Attribute 'densityMin' was not defined");
  if(m_rhomax == HUGE_VAL)
    throw gError("PressureCalculation::setup", "Attribute 'densityMax' was not defined");
  if(m_arraysize_temperature == 0)
    throw gError("PressureCalculation::setup", "Attribute 'arraysize_temperature' was not defined.");
  if(m_arraysize_density == 0)
    throw gError("PressureCalculation::setup", "Attribute 'arraysize_density' was not defined.");

  if(m_Tmin >= m_Tmax)
    throw gError("PressureCalculation::setup", "Attribute 'temperatureMin' (= " + ObjToString(m_Tmin) + ") >= attribute 'temperatureMax' (= " + ObjToString(m_Tmax) + ") not meaningful!");

  if(m_rhomin >= m_rhomax)
    throw gError("PressureCalculation::setup", "Attribute 'densityMin' (= " + ObjToString(m_rhomin) + ") >= attribute 'densityMax' (= " + ObjToString(m_rhomax) + ") not meaningful!");
  
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
	throw gError("PressureCalculation::setup", "Symbol '" + m_symbolName + "' was already created by other module.");
      }
      else
        m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;

      if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName)) {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
          throw gError("PressureCalculation::setup", "Symbol '" + m_temperatureName + "' already exists as a non-scalar.");
        else m_temperatureOffset = Particle::s_tag_format[m_colour].offsetByName(m_temperatureName);
      } else 
          throw gError("PressureCalculation::setup", "Symbol '" + m_temperatureName + "' does not exist but required by this module.");

      if(Particle::s_tag_format[m_colour].attrExists(m_densityName)) {
        if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_densityName).datatype)
          throw gError("PressureCalculation::setup", "Symbol '" + m_densityName + "' already exists as a non-scalar.");
        else m_densityOffset = Particle::s_tag_format[m_colour].offsetByName(m_densityName);
      } else 
          throw gError("PressureCalculation::setup", "Symbol '" + m_densityName + "' does not exist but required by this module.");

      // is it the last cache to be created?
      if(m_colour == M_MANAGER->nColours()-1)
      {
        if(m_phaseUser == 0)
          Particle::registerCache_0(this);
        else if(m_phaseUser == 1) 
          Particle::registerCache(this);
        else // so it is 2
        {
          ParticleCache* pc = new PressureCalculation(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((PressureCalculation*) pc)->stage() == m_stage);
          assert(((PressureCalculation*) pc)->m_colour == m_colour);
          assert(((PressureCalculation*) pc)->m_offset == m_offset);
	  assert(((PressureCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	  assert(((PressureCalculation*) pc)->m_densityOffset == m_densityOffset);
          Particle::registerCache(pc);
          Particle::registerCache_0(this);
        }
      }
      // No? Then make a copy
      else 
      {
        ParticleCache* pc = new PressureCalculation(*this);
        assert(pc->mySymbolName() == m_symbolName);
        assert(((PressureCalculation*) pc)->stage() == m_stage);
        assert(((PressureCalculation*) pc)->m_colour == m_colour);
        assert(((PressureCalculation*) pc)->m_offset == m_offset);
	assert(((PressureCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	assert(((PressureCalculation*) pc)->m_densityOffset == m_densityOffset);

        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);
          
          pc = new PressureCalculation(*this);
          assert(pc->mySymbolName() == m_symbolName);
          assert(((PressureCalculation*) pc)->stage() == m_stage);
          assert(((PressureCalculation*) pc)->m_colour == m_colour);
          assert(((PressureCalculation*) pc)->m_offset == m_offset);
	  assert(((PressureCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
	  assert(((PressureCalculation*) pc)->m_densityOffset == m_densityOffset);
          
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
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_symbolName).datatype)
        throw gError("ParticleCacheDensitySelfContribution::setup", "Symbol '" + m_symbolName + "' was already created by other module.");
    }
    else
      m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false/*persistency*/, m_symbolName).offset;

    if(Particle::s_tag_format[m_colour].attrExists(m_temperatureName)) {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_temperatureName).datatype)
        throw gError("PressureCalculation::setup", "Temperature " + m_temperatureName + " already exists as a non-scalar.");
      else m_temperatureOffset = Particle::s_tag_format[m_colour].offsetByName(m_temperatureName);
    } else 
      throw gError("PressureCalculation::setup", "Symbol '" + m_temperatureName + "' does not exist but required by this module.");

    if(Particle::s_tag_format[m_colour].attrExists(m_densityName)) {
      if(m_datatype != Particle::s_tag_format[m_colour].attrByName(m_densityName).datatype)
        throw gError("PressureCalculation::setup", "Temperature " + m_densityName + " already exists as a non-scalar.");
      else m_densityOffset = Particle::s_tag_format[m_colour].offsetByName(m_densityName);
    } else 
      throw gError("PressureCalculation::setup", "Symbol '" + m_densityName + "' does not exist but required by this module.");

    if(m_phaseUser == 0)
      Particle::registerCache_0(this);
    else if(m_phaseUser == 1)
      Particle::registerCache(this);
    else
    {
      ParticleCache* pc = new PressureCalculation(*this);
      assert(pc->mySymbolName() == m_symbolName);
      assert(((PressureCalculation*) pc)->stage() == m_stage);
      assert(((PressureCalculation*) pc)->m_colour == m_colour);
      assert(((PressureCalculation*) pc)->m_offset == m_offset);
      assert(((PressureCalculation*) pc)->m_temperatureOffset == m_temperatureOffset);
      assert(((PressureCalculation*) pc)->m_densityOffset == m_densityOffset);

      Particle::registerCache(pc);
      Particle::registerCache_0(this);
    } 
  }

  // setup the look-up table
  setupLUT();
}


void PressureCalculation::registerWithParticle()
{
}

