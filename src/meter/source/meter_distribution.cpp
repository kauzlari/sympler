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



#include "meter_distribution.h"

#include "simulation.h"
#include "data_format.h"
#include "manager_cell.h"

using namespace std;

/* Register this Meter with the factory. */
const Meter_Register<MeterDistribution> meter_distribution("MeterDistribution");

#define M_SIMULATION ((Simulation*) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

//---- Constructors/Destructor ----

MeterDistribution::MeterDistribution(Simulation *simulation)
  : Meter(simulation)
{
  init();
}


MeterDistribution::MeterDistribution(Simulation* simulation, const size_t& everyN/*, bool only*/)
  : Meter(simulation, everyN/*, only*/)
{
  init();
}


/*virtual*/ MeterDistribution::~MeterDistribution()
{

}



//---- Methods ----


#define FILL_BINS(var)                                               \
{                                                                    \
  m_n_particles = 0;						     \
  FOR_EACH_PARTICLE_C                                                \
    (phase, m_colour,                                                \
     int bin = (int) (((var - m_min) / m_bin_spacing) + 0.5);        \
     if (bin >= 0 && bin < m_nbins)                                  \
       m_bins->operator[](bin)++;                                    \
     m_n_particles++;						     \
    );                                                               \
} while(0)


/*virtual*/ void MeterDistribution::measureNow(const double& time)
{
  Phase *phase = ((Simulation*) m_parent)->phase();

  if (!m_step) {
    if (!m_bins.isNull()) {
      //      double sum = 0;
      data_sp d = m_format.newData();

      //      for (int i = 0; i < m_nbins; i++)
      //        sum += (*m_bins)[i]*m_bin_spacing;

      for (int i = 0; i < m_nbins; i++)
	(*m_bins)[i] /= m_n_particles*m_bin_spacing*m_avg_over;
	//        (*m_bins)[i] /= sum;

      d->doubleByIndex(0) = time;
      d->vectorDoubleByIndex(1) = m_bin_positions;
      d->vectorDoubleByIndex(2) = m_bins;

      distribute(d);

      /* Don't overwrite, release in case an output module is still holding
         our smart pointer. */
      m_bins.release();
    }

    m_bins.alloc();
    m_bins->resize(m_nbins);
  }


  if (m_offset == -1) {
    if (m_what == "velX") {
      FILL_BINS(__iSLFE->v.x);
    } else if (m_what == "velY") {
      FILL_BINS(__iSLFE->v.y);
    } else if (m_what == "velZ") {
      FILL_BINS(__iSLFE->v.z);
    }
  } else {
    FILL_BINS(__iSLFE->tag.doubleByOffset(m_offset));
  }

  m_step++;
  if (m_step == m_avg_over)
    m_step = 0;
}


void MeterDistribution::setup()
{
  Meter::setup();

  m_colour = M_MANAGER->getColour(m_species);

  m_offset = -1;

  if (m_what != "velX" && m_what != "velY" && m_what != "velZ") {
    if (Particle::s_tag_format[m_colour].attrExists(m_what)) {
      DataFormat::attribute_t attr;

      attr = Particle::s_tag_format[m_colour].attrByName(m_what);

      if (attr.datatype != DataFormat::DOUBLE)
        throw gError
          ("MeterDistribution::setup",
           "Please specify a double degree of freedom.");

      m_offset = attr.offset;
    } else
      throw gError
        ("MeterDistribution::setup",
         "Unknown degree of freedom: '" + m_what + "'. " +
         "Possibilities are: " + Particle::s_tag_format[m_colour].toStr/*ing*/());
  }

  m_size = m_max - m_min;
  m_step = 0;
  m_bin_spacing = m_size / (m_nbins-1);

  m_bin_positions.alloc();
  m_bin_positions->resize(m_nbins);

  for (int i = 0; i < m_nbins; i++)
    m_bin_positions->operator[](i) = m_min + i * m_bin_spacing;
}


void MeterDistribution::init()
{
  m_properties.setClassName("MeterDistribution");

  m_properties.setDescription("Determines the distribution of a particles degree of freedom.");

  STRINGPC
    (species, m_species,
     "Species for measurement.");

  STRINGPC
    (what, m_what,
     "What to measure: Possibilities are 'velX', 'velY', 'velZ' or any other degree of freedom.");

  INTPC
    (nBins, m_nbins, 0, "Number of bins to use for distribution generation.");

  m_properties.addProperty
    ("min",
     PropertyList::DOUBLE,
     &m_min,
     NULL,
     "Minimum distribution value.");

  m_properties.addProperty
    ("max",
     PropertyList::DOUBLE,
     &m_max,
     NULL,
     "Maximum distribution value.");

  INTPC
    (avgOver, m_avg_over, 0, "Average over timesteps.");

  m_nbins = 10;
  m_avg_over = 1;
  m_min = 0;
  m_max = 1;
  m_what = "velX";
  m_species = "UNDEF";

  m_format.addAttribute("time", DataFormat::DOUBLE);
  m_format.addAttribute("value", DataFormat::VECTOR_DOUBLE);
  m_format.addAttribute("probability", DataFormat::VECTOR_DOUBLE);
}
