/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


#include "pca_matrix_inverse.h"
#include "simulation.h"
#include "manager_cell.h"
#include "colour_pair.h"
#include "triplet_calculator.h"
#include "quintet_calculator.h"


const SymbolRegister<PCaMatrixInverse> pca_matrix_inverse("MatrixInverse");

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCaMatrixInverse::PCaMatrixInverse(/*Node*/Simulation* parent)
  : ParticleCache(parent)
{
  m_datatype = DataFormat::TENSOR;
  m_working_mat = gsl_matrix_alloc(SPACE_DIMS, SPACE_DIMS);

  m_inverse.size1 = m_inverse.size2 = m_inverse.tda = SPACE_DIMS;
  m_inverse.owner = 0;
  m_inverse.data = NULL;
  
  m_permutation = gsl_permutation_alloc(3);
    
  init(); 
}

PCaMatrixInverse::~PCaMatrixInverse()
{
  gsl_matrix_free(m_working_mat);
  gsl_permutation_free(m_permutation);
}

void PCaMatrixInverse::init()
{
  m_properties.setClassName("PCaMatrixInverse");
  m_properties.setName("MatrixInverse");

  m_properties.setDescription("Computes the inverse of a given symmetric square matrix.");
      
  STRINGPC
      (symbol, m_symbolName,
       "Symbol name for the inverse.");
  
  STRINGPC
      (tensor, m_tensor_symbol,
       "Name of the tensor for which to compute the inverse.");
  
  m_symbolName = "undefined";
  m_tensor_symbol = "undefined";   
}

void PCaMatrixInverse::setup()
{
  // turns of the default gsl_error handling so that we can handle ourselves
  gsl_set_error_handler_off ();

  if(m_species == "undefined")
    throw gError("PCaMatrixInverse::setup", "Attribute 'species' has value \"undefined\"."); 

  if(m_symbolName == "undefined")
    throw gError("PCaMatrixInverse::setup", "Attribute 'symbol' has value \"undefined\"."); 

  if(m_tensor_symbol == "undefined")
    throw gError("PCaMatrixInverse::setup", "Attribute 'tensor' has value \"undefined\"."); 

  if(m_phaseUser != 0 && m_phaseUser != 1 && m_phaseUser != 2)
    throw gError("PCaMatrixInverse::setup", "Attribute 'stage' has none of the allowed values \"0\", \"1\", \"2\".");
  
  pair<size_t, size_t> tempPair;
  
  // should we create a Cache for the other colours too?
  if(m_species == "ALL")
  {
    for (m_colour = 0; m_colour < M_MANAGER->nColours(); ++m_colour)
    {
      // is the symbol already existing somewhere?
      if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
        throw gError("PCaMatrixInverse::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed in MatrixInverse.");
      // the tensor MUST already exist
      if(!Particle::s_tag_format[m_colour].attrExists(m_tensor_symbol))
        throw gError("PCaMatrixInverse::setup", "Symbol " + m_tensor_symbol + " not found for species '" + M_MANAGER->species(m_colour) + "'.");
    }
    m_colour = 0;
    if(Particle::s_tag_format.size() > 1)
    {
      // loop over colours, except the last one for making copies
      for(m_colour = 0; m_colour < Particle::s_tag_format.size()-1; ++m_colour)
      {
        // FIXME: persistency currently always false -> is this OK?
        m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
        m_tensor_offset = Particle::s_tag_format[m_colour].offsetByName(m_tensor_symbol);
                
        // register a copy
        ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
        assert(((PCaMatrixInverse*) pc)->mySymbolName() == m_symbolName);
        assert(((PCaMatrixInverse*) pc)->stage() == m_stage);
        assert(((PCaMatrixInverse*) pc)->m_colour == m_colour);
        assert(((PCaMatrixInverse*) pc)->m_offset == m_offset);
        
        if(m_phaseUser == 0)
          Particle::registerCache_0(pc);
        else if(m_phaseUser == 1)
          Particle::registerCache(pc);
        else // so it is 2
        {
          Particle::registerCache(pc);
          
          pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
          assert(((PCaMatrixInverse*) pc)->mySymbolName() == m_symbolName);
          assert(((PCaMatrixInverse*) pc)->stage() == m_stage);
          assert(((PCaMatrixInverse*) pc)->m_colour == m_colour);
          assert(((PCaMatrixInverse*) pc)->m_offset == m_offset);
           
          Particle::registerCache_0(pc);
        }
      }
      // now, m_colour = last colour 
    }
  }
  else
  {
    m_colour = M_MANAGER->getColour/*AndAdd*/(m_species);
    // are the symbols already existing?
    if(Particle::s_tag_format[m_colour].attrExists(m_symbolName))
      throw gError("PCaMatrixInverse::setup", "Symbol " + m_symbolName + " is already existing for species '" + M_MANAGER->species(m_colour) + "'. Second definition is not allowed in MatrixInverse.");
    // the tensor MUST already exist
    if(!Particle::s_tag_format[m_colour].attrExists(m_tensor_symbol))
      throw gError("PCaMatrixInverse::setup", "Symbol " + m_tensor_symbol + " not found for species '" + M_MANAGER->species(m_colour) + "'.");
  }
  // next lines are done in any case
  // FIXME: persistency currently always false -> is this OK?
  m_offset = Particle::s_tag_format[m_colour].addAttribute(m_symbolName, m_datatype, false, m_symbolName).offset;
  
//   MSG_DEBUG("PCaMatrixInverse::setup", "m_evecs_offset = " << m_evecs_offset);
  
  m_tensor_offset = Particle::s_tag_format[m_colour].offsetByName(m_tensor_symbol);
  
  if(m_phaseUser == 0) 
    Particle::registerCache_0(this);
  else if(m_phaseUser == 1) 
    Particle::registerCache(this);
  else // so it is 2
  {
    // register a copy
    ParticleCache* pc = copyMySelf()/*new ParticleCacheArbitrary(*this)*/;
    assert(((PCaMatrixInverse*) pc)->mySymbolName() == m_symbolName);
    assert(((PCaMatrixInverse*) pc)->stage() == m_stage);
    assert(((PCaMatrixInverse*) pc)->m_colour == m_colour);
    assert(((PCaMatrixInverse*) pc)->m_offset == m_offset);
        
    Particle::registerCache(pc);
    Particle::registerCache_0(this);
    
  }
}


bool PCaMatrixInverse::findStage()
{
  return Symbol::findStageNewPrelim();
}


bool PCaMatrixInverse::findStage_0()
{
  return Symbol::findStageNewPrelim_0();
}
