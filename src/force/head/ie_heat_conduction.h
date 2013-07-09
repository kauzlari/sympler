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



#ifndef __IE_HEAT_CONDUCTION_H
#define __IE_HEAT_CONDUCTION_H 

#include "f_with_rng.h"
#include "pairdist.h"
#include "manager_cell.h"
#include "integrator_energy.h"

using namespace std;


#define KAPPA 1

//---- IEHeatConduction ----

/*!
 * Mesoscopic heat conduction term from the DPDE model. See
 * J. B. Avalos and A. D. Mackie, Europhys. Lett. 40, 141-146 (1997) and
 * P. Espanol, Europhys. Lett. 40, 631-636 (1997)
 */
class IEHeatConduction : public FWithRng
{
protected:
  /*!
   * Tag offset of the force on the internal energy
   */
  pair<int, int> m_force_offset[FORCE_HIST_SIZE];

  /*!
   * Tag offset of the internal energy
   */
  pair<int, int> m_energy_offset;

  /*!
   * Reciprocal of the cut-off radius
   */
  double m_rcinv;

  /*!
   * = 1/sqrt(dt)
   */
  double m_r_sqrt_dt;

  /*!
   * Mesoscopic noise amplitude for the heat conduction term
   */
  FunctionFixed m_alpha;

  /*!
   * Mesoscopic heat conduction coefficient
   */
  FunctionFixed m_kappa;

  /*!
   * Tag offset of the reciprocal of the pair distance
   */
  size_t m_compute_ri_offset;

  /*!
   * Integrator for the energy degree of freedom. Needed for temperature calculation.
   */
  pair<IntegratorEnergy*, IntegratorEnergy*> m_integrator;

  void init();

  /*!
   * Return the noise amplitude for the pair of particles with
   * internal energies \a energy1 and \a energy2, respectively.
   * @param energy1 Internal energy of the first particle
   * @param energy2 Internal energy of the second particle
   */
  double alpha(double energy1, double energy2) {
    return m_alpha(energy1+energy2);
  }

  /*!
   * Return the dissipation for the pair of particles with
   * internal energies \a energy1 and \a energy2, respectively.
   * @param energy1 Internal energy of the first particle
   * @param energy2 Internal energy of the second particle
   */
  double kappa(double energy1, double energy2) {
    return m_kappa(energy1+energy2);
  }

public:
  /*!
   * Constructor
   * @param simulation Pointer to the simulation object
   */
  IEHeatConduction(Simulation *simulation);

  /*!
   * Destructor
   */
  virtual ~IEHeatConduction();

#ifdef _OPENMP
  virtual void setForceSlots(Integrator* intr, int thread_no);
#endif

  virtual void setup();

  virtual void computeForces(int force_index);

  virtual void computeForces(Particle* part, int force_index);

#ifndef _OPENMP
  virtual void computeForces(Pairdist* pair, int force_index);
#else
  virtual void computeForces(Pairdist* pair, int force_index, int thread_no);

//   virtual void mergeCopies(Particle* p, size_t thread_no, int force_index) {}
#endif
};

#endif
