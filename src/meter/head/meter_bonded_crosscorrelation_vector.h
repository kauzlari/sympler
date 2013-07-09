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



#ifndef __METER_BONDED_CROSSCORRELATION_VECTOR_H
#define __METER_BONDED_CROSSCORRELATION_VECTOR_H

//#include "gen_f.h"
#include "meter.h"
#include "colour_pair.h"

using namespace std;

/* --- MeterVisco --- */

class Simulation;

/*!
 * This meter computes the cross-correlation function for two arbitrary vector-propertise V1(t=0) and V2(t) stored in the tag of two different particles. The outer product of the vectors is computed resulting in a tensor. The two particles belong to a list of bonded pairs specified by the user. 
 *
 * ASSUMPTIONS:
 * - constant timestep length
 * - the order of the particles in the pairs of the connected list is always the same with respect to the species and other important properties the user might think of. For example the connections might all be pointing from south to north with the first particle being in the south and the second one being in the north.
 */
class MeterBondedCrosscorrelationVector: public Meter
{
protected:

  /*!
   * index of the bonded list, this calculator belongs to
   */
  size_t m_listIndex;
    
  /*!
   * name of the bonded list, this calculator belongs to
   */
  string m_listName;

  /*!
   * The two species containing the input data.
   */
  pair<string, string> m_species;

  /*!
   * The two species containing the input data.
   */
  pair<size_t, size_t> m_colours;

  /*!
   * The \a ColourPair where we find the bonded list this 
   * \a MeterBondedCrosscorrelationVector is working on
   */
  ColourPair* m_cp;


  /*!
   * storage for individual autocorrelation functions to be averaged later in \a m_cfAv
   */
  vector<vector<tensor_t> > m_cf;

  /*!
   * storage for the autocorrelation function averaged from \a m_cf
   */
  vector<tensor_t> m_cfAv;
  
  /*!
   * memory offset to the origin (t=0-value) of each currently measured \a m_cf . The origin is stored in the \a Particle tag 
   */
  /*pair<size_t,*/ size_t /*>*/ m_cfOrgOffset; // trying with offset to VECTOR_INT

  /*!
   * Symbols of the two input values for which the CF should be computed
   */
  pair<string, string> m_inputSymbol;

  /*!
   * memory offset to the input values of the two \a m_species for which the CF should be computed 
   */
  pair<size_t, size_t> m_inputOffset;

  /*!
   * current index for each of the acfs stored in \a m_cf
   */
  vector<int> m_indexCf;

  /*!
   * The number of completed contributions from \a m_cf to \a m_cfAv
   */
  size_t m_countCfAv;

  /*!
   * The number of contributions from \a m_cf to \a m_cfAv to be computed
   */
  size_t m_limitCfAv;

  /*!
   * The number of buffers for parallel accumulation of data in \a m_cf
   */
  size_t m_nBuffCf;

  /*!
   * The length of the correlation function
   */
  size_t m_nValCf;


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
  void zeroCf();

  /*!
   * If complete, collect the individual contributions in \a m_acf into \a m_cfAv.
   * If we have \a m_limitCfAv contributions completed, compute and output the 
   * final result.
   */
  void accumCf(/*const double& time*/);


public:

  /*!
   * Constructor
   */
  MeterBondedCrosscorrelationVector(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~MeterBondedCrosscorrelationVector();
	
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

