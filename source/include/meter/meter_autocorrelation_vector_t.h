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



#ifndef __METER_AUTOCORRELATION_VECTOR_T_H
#define __METER_AUTOCORRELATION_VECTOR_T_H

#include "gen_f.h"
#include "meter.h"

using namespace std;

/* --- MeterVisco --- */

class Simulation;

/*!
 * This meter computes the tensorial autocorrelation funtion for two arbitrary vector-properties V1(t=0) and V2(t) stored in the \a Particle 's tag. Note that the autocorrelation of an average vector quantity of the WHOLE SYSTEM CANNOT be computed by this \a Meter !
 *
 * ASSUMPTIONS:
 * - constant timestep length
 *
 */
class MeterAutocorrelationVectorT: public Meter
{
protected:
  /*!
   * Species containing the input data.
   */
  string m_species;

  /*!
   * Integer identifier corresponding to \a m_species
   */
  size_t m_colour;

  /*!
   * storage for individual autocorrelation functions to be averaged later in \a m_acfAv
   */
  vector<vector<tensor_t> > m_acf;

  /*!
   * storage for the autocorrelation function averaged from \a m_acf
   */
  vector<tensor_t> m_acfAv;
  
  /*!
   * memory offset to the origin of each currently measured \a m_acf . The origin is stored in the \a Particle tag 
   */
  /*vector<*/ size_t /*>*/ m_acfOrgOffset; // trying with offset to VECTOR_INT

  /*!
   * Symbol of the input value for which the ACF should be computed
   */
  pair<string, string> m_inputSymbol;

  /*!
   * memory offset to the input values V1(t=0) and V2(t) for which the ACF should be computed 
   */
  pair<size_t, size_t> m_inputOffset;

  /*!
   * current index for each of the acfs stored in \a m_acf
   */
  vector<int> m_indexAcf;

  /*!
   * The number of completed contributions from \a m_acf to \a m_acfAv
   */
  size_t m_countAcfAv;

  /*!
   * The number of contributions from \a m_acf to \a m_acfAv to be computed
   */
  size_t m_limitAcfAv;

  /*!
   * The number of buffers for parallel accumulation of data in \a m_acf
   */
  size_t m_nBuffAcf;

  /*!
   * The length of the autocorrelation function
   */
  size_t m_nValAcf;


  /*!
   * Storage for the time axis of the acf
   */
  vector<double> m_time;

  /*!
   * helper stating whether the time axis of the acf in \a m_time is already complete
   */
bool  m_timeDone;

  /*!
   * Initialise the property list
   */
  void init();

  /*!
   * Reset the measurement
   */
  void zeroAcf();

  /*!
   * If complete, collect the individual contributions in \a m_acf into \a m_acfAv.
   * If we have \a m_limitAcfAv contributions completed, compute and output the 
   * final result.
   */
  void accumAcf(/*const double& time*/);


public:

  /*!
   * Constructor
   */
  MeterAutocorrelationVectorT(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~MeterAutocorrelationVectorT();
	
  /*!
   * The measurement
   */
  virtual void measureNow(const double& time);
		
  /*!
   * Setup this \a Meter
   */
  virtual void setup();
};

#endif

