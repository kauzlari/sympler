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



#ifndef __METER_LIVE_GEOMVIEW_H_
#define __METER_LIVE_GEOMVIEW_H_

#include "config.h"
#include "meter.h"
#include <stdio.h>

using namespace std;

class Simulation;

class MeterLiveGeomview : public Meter
{
protected:

  /*!
   * The name of the geomview executable
   */
  string m_geomexec;

  /*!
   * The pipe to connect to geomview
   */
  FILE *m_pipe;

  /*!
   * Initialize the property list
   */
  void init();
  
  /*!
   * Convert litte endian float to big endian float
   */
  float htonf(float f);

public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  MeterLiveGeomview(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~MeterLiveGeomview();

  /*!
   * Initialize the geomview display
   */
  virtual void setup();

  /*!
   * Draw picture on screen
   */
  virtual void measureNow(const double& time);
};


#endif /*__METER_LIVE_GEOMVIEW_H_*/
