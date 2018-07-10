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



#ifndef __GRID_AVERAGER_H
#define __GRID_AVERAGER_H

#include "meter.h"
#include "grid_meter.h"
#include "data_format.h"

using namespace std;



/*!
 * Subdivide the simulation domain in small cells and calculate average for
 * each of those cells.
 */
class GridAverager: public Meter
{
 protected:
  /*!
   * Stores the number of steps to average over
   */
  size_t m_avg_steps;

  /*!
   * Current step
   */
  size_t m_step;

  /*!
   * Total number of cells
   */
  int m_n_cells;

  /*!
   * Index where the data from the GridMeters starts
   */
  int m_idx_data_start;

  /*!
   * Which cell does a certain particle lie in?
   */
  vector<vector<size_t> > m_locations;

  /*!
   * Number of particles per cell
   */
  vector<vector<size_t> > m_n_particles;

  /*!
   * Number of total particles per cell (independent of color)  
   */
  vector<size_t> m_n_total_particles;

  /*!
   * List of all \a GridMeter s, which determine the quantities within each cell
   */
  list<GridMeter*> m_meters;
    
  /*!
   * Initialize the property list
   */
  void init();

  /*!
   * Calculate averages for this step
   */
  void average(data_sp data);

  /*!
   * Finish this averaging step
   */
  void finishStep(data_sp data);

  /*!
   * Load GridMeters from XML
   */
  virtual Node* instantiateChildWithXML(const string &name, const xmlNode *xmln);

  /*!
   * Write body of the XML file
   */
  virtual void writeMoreBody(ostream &s, int shift);

public:
  /*!
   * Constructor
   * @param simulation Pointer to the parent \a Simulation object
   */
  GridAverager(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~GridAverager();

  /*!
   * Initialize the location arrays and setup child \a GridMeter s
   */
  virtual void setup();

  
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
      // measure if we are in an averaging cycle
      if (m_step) {
        ++m_counter;
        measureNow(time);
      }
      // the case m_step == m_counter is excluded by construction
      else if(m_counter++ == m_measure_every_n)
      {
        m_counter = 1;
        measureNow(time);
      }
    }
  }
		

  
  /*!
   * Return location of the particle with color \a colour and index \a i
   * @param colour Color of the particle
   * @param i Index (slot) of the particle in the \a Phase
   */
  size_t location(size_t colour, size_t i) const {
    return m_locations[colour][i];
  }

  /*!
   * Return number of particles of color \a colour in cell \a i
   * @param colour Color of the particle
   * @param i Index of the cell
   */
  size_t nParticles(size_t colour, int i) const {
    return m_n_particles[colour][i];
  }

  /*!
   * Return total number of particles in cell \a i
   * @param i Index of the cell
   */
  size_t nParticles(int i) const {
    return m_n_total_particles[i];
  }

  /*!
   * Return the volume of cell \a i
   * @param i Index of the cell
   */
  virtual double volume(int i) const = 0;

  /*!
   * Return the number of cells
   */
  int nCells() {
    return m_n_cells;
  }

  /*!
   * Return the format of the output of this meter
   */
  DataFormat &format() {
    return m_format;
  }
};


#endif
