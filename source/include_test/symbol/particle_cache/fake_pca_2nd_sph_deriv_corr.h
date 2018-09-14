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


#ifndef __FAKE_PCA_2ND_DERIV_CORR_H
#define __FAKE_PCA_2ND_DERIV_CORR_H

#include "pca_2nd_sph_deriv_corr.h"
#include "fake_particle_cache.h"

/*!
 * Fake class for testing class \a PCa2ndSPHDerivCorr
 */
class FakePCa2ndSPHDerivCorr: public PCa2ndSPHDerivCorr
{

  protected:


  public:

    /*!
     * Constructor for node hierarchy
     */
    FakePCa2ndSPHDerivCorr(/*Node*/Simulation* parent);

    /*!
     * Destructor
     */
    virtual ~FakePCa2ndSPHDerivCorr();

    /*!
     * Set protected member \a PCa2ndSPHDerivCorr::m_systemMatOffset
     */
    virtual void setSystemMatOffset(size_t offset) {

    	m_systemMatOffset = offset;
    }

    /*!
     * Set protected member \a PCa2ndSPHDerivCorr::m_offset
     */
    virtual void setOutputOffset(size_t offset) {

    	m_offset = offset;
    }

    /*!
     * Set protected member \a PCa2ndSPHDerivCorr::m_pairLoopToDo
     */
    virtual void setPairLoopToDo(bool pairLoopToDo) {

    	m_pairLoopToDo = pairLoopToDo;
    }

    /*!
     * Set protected member \a PCa2ndSPHDerivCorr::m_colour
     */
    virtual void setColour(size_t colour) {

    	m_colour = colour;
    }

};

#endif


