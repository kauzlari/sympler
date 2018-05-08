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


#include "cbl_pair_particle_arbitrary.h"

// implicitly needed but included by colour_pair.h
// #include "simulation.h"
// #include "manager_cell.h"
#include "colour_pair.h"


#define M_SIMULATION  ((Simulation*) m_parent)
#define M_PHASE  M_SIMULATION->phase()
#define M_MANAGER  M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

CblPairParticleArbitrary::CblPairParticleArbitrary
(/*Node*/Simulation* parent)
  : Callable(parent)
{
  init();
}    

void CblPairParticleArbitrary::init()
{
  m_properties.setClassName("CblPairParticleArbitrary");

  m_properties.setDescription
    ("Computes a completely user-defined property for each particle, "
     "which needs pair summation. The definition is done with the "
     "attribute 'expression' and the particleFactors.");
  
  STRINGPC
    (species1, m_species.first,
     "Name for the first species of the pairs, this Symbol is used "
     "for.");
  
  STRINGPC
    (species2, m_species.second,
     "Name for the second species of the pairs, this Symbol is used "
     "for.");
  
  m_species.first = "undefined";
  m_species.second = "undefined";
  
  STRINGPC
    (symbol, m_symbolName,
     "Name of the symbol for the calculated property.");
  
  STRINGPC
    (expression, m_expression,
     "The mathematical expression to be computed for the pairFactor.");
  
  STRINGPC
    (particleFactor_i, m_1stPExpression,
     "The mathematical expression of the additional particle factor "
     "for the first particle.");
  
  STRINGPC
    (particleFactor_j, m_2ndPExpression,
     "The mathematical expression of the additional particle factor "
     "for the second particle.");
  
  INTPC
    (symmetry, m_symmetry, -2,
     "How does the pair expression behave under index interchange? "
     "\"1\" means symmetric behaviour. \"-1\" means antisymmetric "
     "behaviour.");
  
  DOUBLEPC
    (cutoff, m_cutoff, 0,
     "Specifies the range of the pair summation.");
  
  BOOLPC
    (overwrite, m_overwrite,
     "Should an already existing symbol with name 'symbolName' be "
     "overwritten?.");
  
  m_expression = "undefined";
  m_1stPExpression = "undefined";
  m_2ndPExpression = "undefined";
  m_symbolName = "undefined";
  m_symmetry = 1;
  m_cutoff = 0.;
  m_overwrite = false;

}

void CblPairParticleArbitrary::setup()
{
  MSG_DEBUG
    ("CblPairParticleArbitrary::setup", "Class name: " + className() +
     ", name: " + name() + ": START, cutoff = " << m_cutoff);    
    
  if(m_cutoff <= 0)
    throw
      gError
      ("CblPairParticleArbitrary::setup", "Class name: " + className()
       + ", name: " + name() + ": Attribute 'cutoff' has no value "
       "> 0!");
     
  if(m_expression == "undefined")
    throw
      gError
      ("CblPairParticleArbitrary::setup", "Class name: " + className()
       + ", name: " + name() + ": Attribute 'expression' has value "
       "\"undefined\"");
     
  if(m_symmetry != -1 && m_symmetry != 1)
    throw
      gError
      ("CblPairParticleArbitrary::setup", "Class name: " + className()
       + ", name: " + name() + ": Attribute 'symmetry' may only be "
       "\"-1\" (antisymmetry) or \"1\" (symmetry).");

  if(m_species.first == "undefined")
    throw
      gError
      ("CblPairParticleArbitrary::setup for module " + className(),
       "Attribute 'species1' has value \"undefined\"."); 

  if(m_species.second == "undefined")
    throw
      gError
      ("CblPairParticleArbitrary::setup for module " + className(),
       "Attribute 'species2' has value \"undefined\"."); 

  
  m_cp = M_MANAGER ->
    cp(M_MANAGER->getColour(m_species.first),
       M_MANAGER->getColour(m_species.second));
    
  m_cp->setCutoff(m_cutoff);
  m_cp->setNeedPairs(true);

  if(m_overwrite) {
    try {
      m_slots.first = Particle::s_tag_format[m_cp->firstColour()]
	.indexOf(m_symbolName, m_datatype);
      m_slots.first = Particle::s_tag_format[m_cp->firstColour()]
	.offsetByIndex(m_slots.first);
    }
    catch(gError& err) {
      throw
	gError
	("CblPairParticleArbitrary::setup for module " + className(),
	 "Search for symbol for species '"
	 + M_MANAGER->species(m_cp->firstColour()) + " failed. The "
	 "message was " + err.message()); 
    }
    try {
      m_slots.second = Particle::s_tag_format[m_cp->secondColour()]
	.indexOf(m_symbolName, m_datatype);
      m_slots.second = Particle::s_tag_format[m_cp->secondColour()]
	.offsetByIndex(m_slots.second);
    }
    catch(gError& err) {
      throw
	gError
	("CblPairParticleArbitrary::setup for module " + className(),
	 "Search for symbol for species '"
	 + M_MANAGER->species(m_cp->secondColour()) + " failed. The "
	 "message was " + err.message()); 
    }
    
  } // end of if(m_overwrite)
  else {
    
    if(Particle::s_tag_format[m_cp->firstColour()]
       .attrExists(m_symbolName)) 
      throw
	gError
	("CblPairParticleArbitrary::setup for module " + className(),
	 "Symbol " + m_symbolName + " is already existing for "
	 "species '" + M_MANAGER->species(m_cp->firstColour()) + "'. "
	 "Second definition is not allowed for overwrite = \"no\"");
    
    if(Particle::s_tag_format[m_cp->secondColour()]
       .attrExists(m_symbolName))
      throw
	gError
	("CblPairParticleArbitrary::setup for module " + className(),
	 "Symbol " + m_symbolName + " is already existing for "
	 "species '" + M_MANAGER->species(m_cp->secondColour()) + "'. "
	 "Second definition is not allowed for overwrite = \"no\"");        
    
    // see CONVENTION5 for rule about persistencies
    m_slots.first = Particle::s_tag_format[m_cp->firstColour()]
      .addAttribute(m_symbolName, m_datatype, /*persistency*/false,
		    m_symbolName).offset;
    
    if(m_cp->firstColour() != m_cp->secondColour()) {
      // see CONVENTION5 for rule about persistencies
      m_slots.second = Particle::s_tag_format[m_cp->secondColour()]
	.addAttribute(m_symbolName, m_datatype, /*persistency*/false,
		      m_symbolName).offset;
    }
    else m_slots.second = m_slots.first;
    
  } // end of else of if(m_overwrite)
  
  m_function.setExpression(m_expression);
  m_function.setColourPair(m_cp);
  m_1stparticleFactor.setExpression(m_1stPExpression);
  m_1stparticleFactor.setColourPair(m_cp);
  m_2ndparticleFactor.setExpression(m_2ndPExpression);
  m_2ndparticleFactor.setColourPair(m_cp);
}

