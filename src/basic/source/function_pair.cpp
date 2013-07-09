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



#include "function_pair.h"

#include "phase.h"
#include "fp_scalar.h"
#include "simulation.h"
#include "colour_pair.h"
#include "manager_cell.h"
#include "val_calculator_rho.h"
#include "val_calculator_kernel.h"
#include "weighting_function.h"

#include "pca_density_0oc.h"
#include "pca_volume_self_contribution.h"

// #include "val_calculator_shear.h"

#include "fp_vector.h"
#include "fp_tensor.h"


#define M_MANAGER m_cp->manager()
#define M_PHASE M_MANAGER->phase()
#define M_SIMULATION ((Simulation*) M_PHASE->parent())


/* FunctionPair */

FunctionPair::FunctionPair()
{
//   MSG_DEBUG("FunctionPair::FunctionPair()", "CALLED");
}


FunctionPair::~FunctionPair()
{
}



class FunctionPairCallback: public UnknownSymbolCallback
{
protected:
  ColourPair *m_cp;

public:
  FunctionPairCallback(ColourPair *cp): m_cp(cp)
  { 
  }

  virtual TypedValue *operator()(const string &e) {
    return NULL;
    
  }
};


void FunctionPair::compile()
{

  FunctionArbitrary::compile();

  assert(m_cp);

  /*
   * Particle-index
   */
  addInt("i", "first_particle", offsetof(Particle, mySlot));
  addInt("j", "second_particle", offsetof(Particle, mySlot));

  /*
   * Relative distance vector
   */
  addPoint(/*"[rij]"*/"rij", "dist", offsetof(dist_t, cartesian));
  /*
   * Absolute relative distance
   */
  addDouble("rij", "dist", offsetof(dist_t, abs));
  /*
   * Position
   */
  addPoint("ri"/*"[ri]"*/, "first_particle", offsetof(Particle, r));
  addPoint(/*"[rj]"*/"rj", "second_particle", offsetof(Particle, r));
  
  /*
   * Velocity
   */
  addPoint("vi"/*"[vi]"*/, "first_particle", offsetof(Particle, v));
  addPoint(/*"[vj]"*/"vj", "second_particle", offsetof(Particle, v));

  /*
   * Add all that's found in the pair tag
   */
  addAllFromDataFormat("pair_tag", &m_cp->tagFormat(), "ij");  

  /*
   * Add all that's found in the particle tag
   */
  addAllFromDataFormat
      ("first_tag", 
       Particle::s_tag_format[/*m_colour*/m_cp->firstColour()].format(), "i");
  addAllFromDataFormat
      ("second_tag", 
       Particle::s_tag_format[/*m_colour*/m_cp->secondColour()].format(), "j");

  m_parser.setCallback(new FunctionPairCallback(m_cp));
  m_parser.parse(m_expression);

  m_compiler.setHeader
    ("void *result, void *dist, void *pair_tag, "
     "void *first_particle, void *second_particle, void *first_tag, void *second_tag");
  m_compiler.setParserAndCompile(&m_parser);
}
