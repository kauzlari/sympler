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



// should be usable for all dimensionalities

#ifndef __METER_H
#define __METER_H

#include <unistd.h>

#include <list>
#include <string>
#include <fstream>
#include <iostream>	// for void printScreen()

#include "general.h"
#include "smart_enum.h"
#include "data_format.h"
#include "postprocessor.h"
#include "node_many_children.h"

using namespace std;

//---- newclass -----------------------------------------------------------



template<typename T>
class math_matrix_t {
 protected:
  int m_n, m_m;
  T *m_data;
	
 public:
  math_matrix_t(int n = SPACE_DIMS, int m = SPACE_DIMS) : m_n(n), m_m(m) {
    m_data = new T[n * m];
  }
    ~math_matrix_t() {
      delete m_data;
    }

    T &operator()(int y, int x) {
      return m_data[y * m_m + x];
    }
    const T &operator()(int y, int x) const {
      return m_data[y * m_m + x];
    }

    point_t operator*(const point_t &p) {
      point_t r;
	
      for (int j = 0; j < m_n; j++) {
	r[j] = 0;
		
	for (int i = 0; i < m_m; i++) {
	  r[j] += operator()(j, i)*p[i];
	}
      }
	
      return r;
    }

    math_matrix_t<T> operator*(const math_matrix_t<T> &a) {
      math_matrix_t<T> r(m_n, a.m_m);
	
      for (int j = 0; j < m_n; j++) {
	for (int i = 0; i < a.m_m; i++) {
	  for (int k = 0; k < m_m; k++) {
	    r(j, i) += operator()(j, k)*a(k, i);
	  }
	}
      }
	
      return r;
    }
};


struct disp_point_t {
  point_t r;
  double p;
};


typedef math_matrix_t<double> matrix_t;


class Simulation;

/*!
 * This is the base class to all \a Meter s, which are modules that extract
 * measurable quantities (usually average) from the simulation.
 */ 
class Meter: public NodeManyChildren
{
protected:
  /*!
   * Number of timesteps between measurements.
   */
  size_t m_measure_every_n;

  /*!
   * Number of steps since the last measurement
   */
  size_t m_counter;

  /*!
   * Start measurement from this step on
   */
  size_t m_from_step_on;

  /*!
   * Start measurement from this time on = m_from_step_on * dt
   */
  double m_from_time_on;

  /*!
   * The format description for this meter
   */
  DataFormat m_format;

  /*!
   * Should be used to initialize the DataFormat class.
   * Initialize the property list
   */
  void init();

  /*!
   * Loads postprocessors
   */
  virtual Node* instantiateChild(const string &name);

  /*!
   * Distribute the measured data to the postprocessors
   */
  inline void distribute(data_sp data) {
    for (list<Node*>::iterator i = m_children.begin(); i != m_children.end(); i++)
      ((Postprocessor*) *i)->push(data);
  }

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  Meter(Simulation *simulation);

  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   * @param everyN Measurement interval
   */
  Meter(Simulation *simulation, size_t everyN/*, bool only*/);

  /*!
   * Destructor
   */
  virtual ~Meter();
	
  /*!
   * This one is being called right before the simulation is about to start.
   */
  virtual void aboutToStart();

  /*!
   * Setup helper variables
   */
  virtual void setup();

  /*!
   * This one is being called at the end of the simulation to write everything to
   * disk, network, etc.
   */
  virtual void flush();

  /*!
   * Perfom the actual measurement. This is called by the \a Controller directly before
   * the integration procedure. Thus, the first measurement will always be the
   * initial condition.
   *
   * \a measure keeps track of the current time step and calls \a measureNow whenever
   * the measurement interval is reached.
   */
  virtual void measure(const double& time) {
    if(time >= m_from_time_on) {
      if (m_counter++ == m_measure_every_n) {
	m_counter = 1;
	measureNow(time);
      }
    }
  }
		
  /*!
   * Perform the measurement
   */
  virtual void measureNow(const double& time) = 0;
};

//---- Factories ----

class Meter_Factory: public SmartEnum<Meter_Factory>
{
 public:
  virtual Meter *instantiate(Simulation *simulation) const = 0;

 protected:
  Meter_Factory(const string &name)
    : SmartEnum<Meter_Factory>(name) { }
};


template <class T>
class Meter_Register: public Meter_Factory
{
 public:
  Meter_Register(const string &name)
    : Meter_Factory(name) { }

    virtual Meter *instantiate(Simulation *simulation) const;
};



//---- Inline functions ----

template <class T>
inline Meter *Meter_Register<T>::instantiate(Simulation *simulation) const
{
  return new T(simulation);
}



#endif  
