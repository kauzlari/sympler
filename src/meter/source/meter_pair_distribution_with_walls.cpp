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



#include "meter_pair_distribution_with_walls.h"

#include "simulation.h"
#include "data_format.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterPairDistributionWithWalls> meter_pair_distribution("MeterPairDistributionWithWalls");

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

MeterPairDistributionWithWalls::MeterPairDistributionWithWalls(Simulation *simulation)
  : Meter(simulation)
{
  init();
}


MeterPairDistributionWithWalls::MeterPairDistributionWithWalls(Simulation* simulation, const size_t& everyN)
  : Meter(simulation, everyN)
{
  init();
}


MeterPairDistributionWithWalls::~MeterPairDistributionWithWalls()
{

}



//---- Methods ----


void MeterPairDistributionWithWalls::measureNow(const double& time)
{
  if (!m_step) {
    if (!m_bins.isNull()) {
      data_sp d = m_format.newData();

      double volume, N, factor;

      volume = M_PHASE->cuboidVolume();
      N = M_PHASE->particles(m_colour).size();

      factor = volume/(N*N);

      for (int i = 0; i < m_nbins; i++) {
// 	double rn = i*m_bin_size;
// 	double rnn = rn+m_bin_size;

	(*m_bins)[i] *= factor/m_avg_over;
      }

      d->doubleByIndex(0) = time;
      d->vectorDoubleByIndex(1) = m_bins;

      distribute(d);

      m_bins.release();
    }

    m_bins.alloc();
    m_bins->resize(m_nbins);
  }

  FOR_EACH_FREE_PAIR
    (m_cp,
     if (pair->abs() < m_max_radius) {
       /* Now we need to calculate the volume of the two slices */
       double V1 = calcVolume(pair->firstPart()->r, pair->abs());
       double V2 = calcVolume(pair->secondPart()->r, pair->abs());

       (*m_bins)[(size_t) (pair->abs()*m_r_bin_size)]+= 1/V1 + 1/V2;
     }
     );

  m_step++;
  if (m_step == m_avg_over)
    m_step = 0;
}


double MeterPairDistributionWithWalls::calcVolume(const point_t &pos, double radius)
{
  size_t bin;
  double rinner, router;
  double dist[2];
  double volume;

  bin = (size_t) (radius*m_r_bin_size);

  rinner = bin * m_bin_size;
  router = (bin+1) * m_bin_size;

  volume = 4*M_PI/3*(router*router*router-rinner*rinner*rinner);

  dist[0] = pos[m_wall_dir] - m_left_wall;
  dist[1] = m_right_wall - pos[m_wall_dir];

  for (int i = 0; i < 2; ++i) {
    if (dist[i] < rinner) {
      /* Completely truncate */
      double h1 = router - dist[i]; /* Outer h */
      double h2 = rinner - dist[i]; /* Inner h */
    
      volume -= M_PI/3 * ( h1*h1 * (3*router - h1) - h2*h2 * (3*rinner - h2) );
    } else if (dist[i] < router) {
      /* Only truncate a small slice */
      double h = router - dist[i];
      
      volume -= M_PI/3 *  h*h * (3*router - h);
    }
  }

  return volume;
}


void MeterPairDistributionWithWalls::setup()
{
  Meter::setup();

  m_colour = M_MANAGER->getColour(m_species);

  m_cp = M_MANAGER->cp(m_colour, m_colour);

  m_cp->setCutoff(m_max_radius);

  m_bin_size = m_max_radius / m_nbins;
  m_r_bin_size = 1/m_bin_size;
  m_step = 0;
}


void MeterPairDistributionWithWalls::init()
{
  m_properties.setClassName("MeterPairDistributionWithWalls");

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
    (wallDir, m_wall_dir, -1,
     "Direction in which to find the wall: 0 = x, 1 = y, 2 = z.");

  m_properties.addProperty
    ("leftWall", PropertyList::DOUBLE, &m_left_wall, NULL,
     "Position of the wall to the left.");

  m_properties.addProperty
    ("rightWall", PropertyList::DOUBLE, &m_right_wall, NULL,
     "Position of the wall to the right.");

  INTPC
    (avgOver, m_avg_over, 0, "Average over timesteps.");

  m_wall_dir = 0;

  m_left_wall = -10;
  m_right_wall = 10;

  m_nbins = 10;
  m_avg_over = 1;
  m_max_radius = 1;
  m_species = "UNDEF";

  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute("value", DataFormat::VECTOR_DOUBLE);
}
