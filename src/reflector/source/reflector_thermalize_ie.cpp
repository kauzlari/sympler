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



#include "reflector_thermalize_ie.h"

#include "phase.h"
#include "boundary.h"
#include "particle.h"
#include "simulation.h"
#include "manager_cell.h"

/* Register this Reflector with the factory. */
const Reflector_Register<ReflectorThermalizeInternalEnergy> reflector_thermalize_internal_energy("ReflectorThermalizeInternalEnergy");


//---- Constructors/Destructor ----

ReflectorThermalizeInternalEnergy::ReflectorThermalizeInternalEnergy(/*Wall *wall*/ Boundary* boundary): ReflectorWithRng(boundary)
{
  init();
}


ReflectorThermalizeInternalEnergy::~ReflectorThermalizeInternalEnergy()
{
}



//---- Methods ----

void ReflectorThermalizeInternalEnergy::init()
{
  m_properties.setClassName
    ("ReflectorThermalizeInternalEnergy");

  m_properties.setDescription
    ("Reflects particles back with a random angle of 90 degrees on "
     "average. The velocity magnitude is changed to a new one, drawn from a Maxwell-Boltzmann "
     "distribution, which corresponds to the temperature defined by 'temperature'. The internal "
     "energy of the particle is set exactly to the value defined by 'internalEnergy'."
     );
    
  FUNCTIONFIXEDPC
    (temperature, m_temperature, 
     "Gives the temperature in dependence on the x, y and z coordinate of where the particle hit.");

  FUNCTIONFIXEDPC
    (internalEnergy, m_internal_energy,
     "Gives the internal energy in dependence on the x, y and z coordinate of where the particle hit.");

  m_temperature.addVariables("x", "y", "z");
  m_temperature_sqrt.addVariables("x", "y", "z");
  m_internal_energy.addVariables("x", "y", "z");
}


void ReflectorThermalizeInternalEnergy::setup()
{
  ReflectorWithRng::setup();

  if (m_temperature.isNull() || m_internal_energy.isNull())
    throw gError
      ("ReflectorThermalizeInternalEnergy::setup",
       "Please provide the temperature and internal energy functions.");

  m_temperature_sqrt.setExpression("sqrt(" + m_temperature.expression() + ")");

  m_ie_offset.resize(Particle::s_tag_format.size());
  for (size_t i = 0; i < m_ie_offset.size(); i++)
    m_ie_offset[i] = Particle::s_tag_format[i].attrByName("internal_energy").offset;
}
