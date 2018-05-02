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



#ifndef __RANDOM_H
#define __RANDOM_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_randist.h>

#define RNG_DEFAULT_SEED 1107

/* Random number generalization. */

class RandomNumberGenerator
{
 protected:
    gsl_rng *m_rng;

 public:
    RandomNumberGenerator();
    virtual ~RandomNumberGenerator();

    /* Generate 0 or 1 with equal probability */
    virtual int binary();

    /* Generate random number from 0 to 1, uniformly distributed */
    virtual double uniform();

    /* Generate normally distributed random number. */
    virtual double normal(double sigma);

    /* Rayleigh distribution. */
    virtual double rayleigh(double sigma);

    virtual void setSeed(size_t seed)
    {
      gsl_rng_set(m_rng, seed);
    }
};


/* The global random number generator. */

//extern RandomNumberGenerator g_rng;


#endif
