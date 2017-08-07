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



#ifndef __PC_FILE_H
#define __PC_FILE_H

#include <fstream>
#include "pc_free_pcalc.h"



//---- Classes ----

class ParticleCreatorFile: public ParticleCreatorFreePCalc
{
protected:
  string m_filename;
  /* list of frozen particles */
  	map<int, ParticleList> m_particles_frozen;
  //is particle inside geometry or not? default yes      
  bool m_particlesInside ; 
  /*!
   * Takes the cutoff from SImulation before the "skin" from the Verlet List is added
   */
  double myCutoff;
  void init();

  static string readNext(ifstream &pos);
    
public:
  ParticleCreatorFile(Boundary *boundary);

  virtual void createParticles();
  virtual void flushParticles();
   /*!
  * In principle for adjusting the box size proposed by the boundary. But in fact, 
  * this ParticleCreator adjusts the box size only indirectly via the second and 
  * third argument
  * @param size the box size proposed by the \a Boundary , this 
  * \a ParticleCreator belongs to 
   * @param frameRCfront Should the 'lower end' of the simulation box be extended 
   * in order to put wall \a Particle s there?
   * @param frameRCend Should the 'upper end' of the simulation box be extended 
   * in order to put wall \a Particle s there?
  */
  virtual void adjustBoxSize(point_t &size, bool_point_t& frameRCfront, bool_point_t& frameRCend);
  void readParticle(Particle &p, ifstream &pos, string &freeOrFrozen );
  virtual void setup();
};

#endif
