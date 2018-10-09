/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
 * David Kauzlaric <david.kauzlaric@imtek.uni-freiburg.de>,
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


#ifndef __PCA_2ND_SPH_DERIV_CORR_TEST_H
#define __PCA_2ND_SPH_DERIV_CORR_TEST_H

#include "fake_pca_2nd_sph_deriv_corr.h"
#include "particle_cacheTest.h"

using namespace std;


class PCa2ndSPHDerivCorrTest : public ParticleCacheTest
{
  CPPUNIT_TEST_SUITE (PCa2ndSPHDerivCorrTest);

  CPPUNIT_TEST (initTest);
  CPPUNIT_TEST (setupTest);
  CPPUNIT_TEST (computeCacheForTest);

  CPPUNIT_TEST_SUITE_END ();
  
 public:

  /*!
   * Initialize objects
   */
  void setUp (void);

  
 protected:
  
  /*!
   * Test \a SymbolFParticleScalar::computeCacheFor(Particle* p)
   */
  void computeCacheForTest (void);

  /*!
   * Test \a SymbolFParticleScalar::init(Particle* p)
   */
  virtual void initTest (void); 

  /*!
   * Test \a PCa2ndSPHDerivCorr::setup()
   */
  void setupTest (void);

  /*!
   * Helper to create and assign to \a m_particleCache the right
   * \a ParticleCache to be tested
   */
  void setUpParticleCache (void);

  /*!
   * Helper for testing correctness of className in \a PropertyList of
   * \a PCa2ndSPHDerivCorr
   */
  virtual void testPropClassName (void);

  /*!
   * Helper for testing correct initialisation of "species" attribute of
   * \a PCa2ndSPHDerivCorr
   */
  virtual void testSpeciesAttr (void);

  
 private:
  
};

#endif
