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



#include "meter_relative_velocity.h"

#include "phase.h"
#include "simulation.h"
#include "data_format.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterRelativeVelocity> meter_relative_velocity("MeterRelativeVelocity");


#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()


//---- Constructors/Destructor

MeterRelativeVelocity::MeterRelativeVelocity(Simulation *simulation)
  : Meter(simulation)
{
  init();
}


/*virtual*/ MeterRelativeVelocity::~MeterRelativeVelocity()
{
}



//---- Methods ----

/*virtual*/ void MeterRelativeVelocity::measureNow(const double& time)
{
//   Phase *phase = M_PHASE;

  double v = 0;
  double vv = 0;
  double vvv = 0;
  double vvvv = 0;
  //  point_t v = {{{ 0, 0, 0 }}};
  size_t n = 0;

  FOR_EACH_PAIR
    (m_cp,
     double rv = pair->cartesian()*(pair->firstPart()->v - pair->secondPart()->v)/pair->abs();

     v += rv;
     vv += rv*rv;
     vvv += rv*rv*rv;
     vvvv += rv*rv*rv*rv;

     n++;
     );

  v /= n;
  vv /= n;
  vvv /= n;
  vvvv /= n;

  /* Calculate moments 1 to 4 */

  data_sp data = m_format.newData();
  data->doubleByIndex(0) = time;
  data->doubleByIndex(1) = v;
  data->doubleByIndex(2) = vv-v*v;
  data->doubleByIndex(3) = vvv-3*v*vv+2*v*v*v;
  data->doubleByIndex(4) = vvvv-4*v*vvv-3*vv*vv+12*v*v*vv-6*v*v*v*v;

  distribute(data);
}


void MeterRelativeVelocity::init()
{
  m_properties.setClassName("MeterRelativeVelocity");

  m_properties.setDescription(
    "Measures the variance of the distribution of relative velocities. "
    "Output columns are 'time' and 'variance'."
  );

  STRINGPC
    (species1, m_species.first,
     "First species.");

  STRINGPC
    (species2, m_species.second,
     "Second species.");

  m_species.first = "UNDEF";
  m_species.second = "UNDEF";

  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute("m1", DataFormat::DOUBLE);
  m_format.addAttribute("m2", DataFormat::DOUBLE);
  m_format.addAttribute("m3", DataFormat::DOUBLE);
  m_format.addAttribute("m4", DataFormat::DOUBLE);
}


void MeterRelativeVelocity::setup()
{
  Meter::setup();

  m_cp = M_MANAGER->cp(M_MANAGER->getColour(m_species.first), M_MANAGER->getColour(m_species.second)/*m_species*/);


}
