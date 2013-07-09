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



/*
  This file contains global constants, like string identifier,
  that are valid for the whole project.
*/

#ifndef __CONSTS_H
#define __CONSTS_H

/* Always align to 2^6 = 64 bit boundaries */
#define DATA_ALIGNMENT 3

#define STR_DELIMITER string("_")
#define STR_FORCE string("force")
#define STR_SHEAR_TENSOR string("shear_rate_tensor")

// #define N_CALCULATOR_STAGES 2


// deprecated: the maximum stage for ValCalculators
// #define VC_MAX_STAGE 1

// the maximum stage for ParticleCaches
// FIXME: deprecated. Now we use Particle::s_maxStage
// remove it if the new way seems to work 
#define PCA_MAX_STAGE 2



#endif
