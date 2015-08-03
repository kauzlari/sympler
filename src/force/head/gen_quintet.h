/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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



#ifndef _G_QUINTET_H
#define _G_QUINTET_H

#include "gen_f.h"
#include <map>
#include <list>
#include <utility>

#include "misc.h"
#include "list"
#include "integrator_position.h"

using namespace std;

class Simulation;

/*!
 * Base class for all triplet forces.
 * It also handles
 * the creation of the triplet list for which the force is relevant.
 */
class GenQuintet : public GenF
{
 protected:
  
  /*!
   * Inititialize the property list.
   */
  void init();
  
  /*! The \a quintetList this \a GenQuintet is acting on. 
   *This pointer must be assigned during setup by requesting 
   *the correct list from the \a Phase
   */
  quintetList* m_QuintetList;
  
 public:
  /*!
   * Constructor
   */
  GenQuintet();	
  GenQuintet(Simulation *simulation);		
  /*!
   * Destructor
   */
  virtual ~GenQuintet();
  
  virtual void setup();

  virtual void setupAfterParticleCreation();

};

#endif







