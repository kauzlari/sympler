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



#include "meter_pair_distribution.h"

#include "simulation.h"
#include "data_format.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterPairDistribution> meter_pair_distribution("MeterPairDistribution");

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

MeterPairDistribution::MeterPairDistribution(Simulation *simulation)
  : Meter(simulation)
{
  init();
}


MeterPairDistribution::MeterPairDistribution(Simulation* simulation, const size_t& everyN)
  : Meter(simulation, everyN)
{
  init();
}


MeterPairDistribution::~MeterPairDistribution()
{

}



//---- Methods ----


void MeterPairDistribution::measureNow(const double& time)
{
  if (!m_step) {
    if (!m_bins.isNull()) {
      data_sp d = m_format.newData();

      double volume, N, factor;

      volume = M_PHASE->cuboidVolume();
      N = M_PHASE->particles(m_colour).size();

      //      factor = volume*m_r_bin_size/(2*M_PI*N*N*m_avg_over);
      factor = 3*volume/(2*M_PI*N*N);

      for (int i = 0; i < m_nbins; i++) {
	double rn = i*m_bin_size;
	double rnn = rn+m_bin_size;

	double f = factor/(rnn*rnn*rnn-rn*rn*rn);

	(*m_bins)[i] *= f/m_avg_over;
	(*m_variances)[i]
	  = sqrt((*m_variances)[i]*f*f/(m_avg_over*(m_avg_over-1))-(*m_bins)[i]*(*m_bins)[i]);

	//        (*m_bins)[i] *= factor/(rn*rn);
      }

      d->doubleByIndex(0) = time;
      d->vectorDoubleByIndex(1) = m_bins;
      d->vectorDoubleByIndex(2) = m_variances;

      distribute(d);

      m_bins.release();
      m_variances.release();
    }

    m_bins.alloc();
    m_bins->resize(m_nbins);

    m_variances.alloc();
    m_variances->resize(m_nbins);
  }

  FOR_EACH_FREE_PAIR
    (m_cp,
     if (pair->abs() < m_max_radius) {
       (*m_bins)[(size_t) (pair->abs()*m_r_bin_size)]++;
     }
     );

  for (int i = 0; i < m_nbins; i++) {
    (*m_variances)[i] = (*m_bins)[i]*(*m_bins)[i];
  }

  m_step++;
  if (m_step == m_avg_over)
    m_step = 0;
}


void MeterPairDistribution::setup()
{
  Meter::setup();

  m_colour = M_MANAGER->getColour(m_species);

  m_cp = M_MANAGER->cp(m_colour, m_colour);

  m_cp->setCutoff(m_max_radius);

  m_bin_size = m_max_radius / m_nbins;
  m_r_bin_size = 1/m_bin_size;
  m_step = 0;
}


void MeterPairDistribution::init()
{
  m_properties.setClassName("MeterPairDistribution");

  m_properties.setDescription
    ("Determines the pair distribution function.");

  STRINGPC
    (species, m_species,
     "Species for measurement.");

  INTPC
    (nBins, m_nbins, 0,
     "Number of bins to use for distribution generation.");

  DOUBLEPC
    (maxRadius, m_max_radius, 0,
     "Maximum radius.");

  INTPC
    (avgOver, m_avg_over, 0, "Average over timesteps.");

  m_nbins = 10;
  m_avg_over = 1;
  m_max_radius = 1;
  m_species = "UNDEF";

  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute("value", DataFormat::VECTOR_DOUBLE);
  m_format.addAttribute("error", DataFormat::VECTOR_DOUBLE);
}
