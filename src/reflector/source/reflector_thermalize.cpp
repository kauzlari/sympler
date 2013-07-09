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



#include "reflector_thermalize.h"

/* Register this Reflector with the factory. */
const Reflector_Register<ReflectorThermalize> reflector_thermalize("ReflectorThermalize");


//gsl_rng *ReflectorThermalize::s_rng = NULL;
//Normal ReflectorThermalize::s_rng_normal;
double c_rt_disp_eps = 1e-10;

//---- Constructors/Destructor ----

ReflectorThermalize::ReflectorThermalize(/*Wall *wall*/ Boundary* boundary): ReflectorWithRng(boundary)
{
  init();
}


ReflectorThermalize::ReflectorThermalize(/*Wall *wall*/ Boundary* boundary, double temperature)
  : ReflectorWithRng(boundary), m_temperature(temperature), m_temperature_sqrt(sqrt(temperature))
{
  init();
}


ReflectorThermalize::~ReflectorThermalize()
{
}



//---- Methods ----

void ReflectorThermalize::init()
{
  m_properties.setClassName("ReflectorThermalize");
  m_properties.setDescription(
    "Reflects particles back with a random angle of 90 degrees on "
    "average. The velocity magnitude is changed to a new one, drawn from a Maxwell-Boltzmann "
    "distribution, which corresponds to the temperature defined by 'temperature'."
  );
    
  DOUBLEPC
    (temperature, m_temperature, 0,
     "Temperature the particles velocity is set to.");

//   m_properties.addProperty
//     ("biasVelocity", PropertyList::DOUBLE, &m_bias_velocity, NULL,
//      "The magnitude of the velocity with which the velocity of the reflected particle will be biased.");

  FUNCTIONFIXEDPC(biasVelocity, m_bias_velocity, 
                  "The magnitude of the velocity with which the velocity of the reflected "
                      "particle will be biased. You can include time dependency with the variable 't' and spatial dependency with the variables 'x', 'y' and 'z'.");

  POINTPC
      (biasDir, m_bias_orientation,
     "Vector used for direction of velocity biasing. How this vector is interpreted depends on the value of the attribute 'rotationalBias'. The magnitude is intended to be defined with the attribute 'biasVelocity'. Non-unit vectors given here, will be normalised.");
  
  BOOLPC
      (rotationalBias, m_rotational_bias,
       "This attribute decides, how the bias velocity is applied. \nIf it is set to 'no', the bias velocity is simply given by biasVelocity * biasVector.\n If it is set to 'yes', the bias velocity depends on the surface orientation according to biasVelocity * (surfaceNormal x secondNormal), where 'surfaceNormal' is the normal to the surface of the wall segment and 'x' is the cross-product."
      );
  
  BOOLPC
      (oneHit, m_oneHit,
        "If this attribute is set to 'true', the tracking of the particles trajectory will be aborted after the first hit of the particle at a wall"
      );
    
  m_temperature = 1;
	m_temperature_sqrt = 1;
//   m_bias_velocity = 0;
  m_bias_orientation.x = m_bias_orientation.y = 0;
  m_bias_orientation.z = 1;
  m_rotational_bias = false;
  m_oneHit = false;

  m_bias_velocity.addVariable("t");
  m_bias_velocity.addVariable("x");
  m_bias_velocity.addVariable("y");
  m_bias_velocity.addVariable("z");
  m_bias_velocity.setExpression("0");

}


void ReflectorThermalize::setup()
{
  ReflectorWithRng::setup();

  m_temperature_sqrt = sqrt(m_temperature);

  m_bias_orientation /= m_bias_orientation.abs();
}
