/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2018, 
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


#include "pca_iapws-if97_1var.h"
#include "simulation.h"
#include "manager_cell.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCacheIAPWSIF97OneVar::PCacheIAPWSIF97OneVar
  (size_t colour, size_t offset, string symbolName)
    : ParticleCache(colour, offset, symbolName), m_LUT(NULL) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  // create space for one pointer to one thermodynamic input variable
  m_inputVarPtrs.resize(1);
}

PCacheIAPWSIF97OneVar::PCacheIAPWSIF97OneVar
    (/*Node*/Simulation* parent)
      : ParticleCache(parent), m_LUT(NULL) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  // create space for one pointer to one thermodynamic input variable
  m_inputVarPtrs.resize(1);

  init();
}


PCacheIAPWSIF97OneVar::~PCacheIAPWSIF97OneVar()
{
  if(m_LUT) {
    delete [] m_LUT;
  }
}

void PCacheIAPWSIF97OneVar::setupLUT() {

  // auxiliary variables
  double var1 = m_var1Min;
  // This vector<double*> should have been resized to 1 in constructor
  m_inputVarPtrs[0] = &(var1);  
  size_t slot;
  
  checkConstraints();
  
  // Step size depending on the size of the array and the range of the input value.
  m_var1StepSize = (m_var1Max - m_var1Min)/(m_arraySizeVar1-1);

  // Allocation of the 2D Array, depending on the given sizes.
  m_LUT = new double[m_arraySizeVar1];

  // Calculate results at support points and store them into the LUT
  for(size_t i = 0; i < m_arraySizeVar1; ++i) {

    freesteamCalculationForState(m_LUT[i]);      
      
    var1 += m_var1StepSize; 
  }
}

void PCacheIAPWSIF97OneVar::calculateResult(double& result) const {

  // assuming m_inputVars[0] has been set correctly beforehand by the
  // caller of this function
  double& inputVar1 = *(m_inputVarPtrs[0]);
  
  // out of bounds?
  if((inputVar1 < m_var1Min) || (inputVar1 > m_var1Max))    
    throw gError("PCacheIAPWSIF97OneVar::calculateResult for module " + className(), m_var1Name + " out of bounds: " + m_var1Name + "=" + ObjToString(inputVar1) + ", admissible [" + m_var1Name + "Min," + m_var1Name + "Max] = [" + ObjToString(m_var1Min) + "," + ObjToString(m_var1Max) + "].");

  // Calculation of the surrounding sampling points in 1D LUT
  size_t var1Slot0 = (floor((inputVar1-m_var1Min)/m_var1StepSize));

  // Normalisation of the input value
  double var1Normalised =
    ((inputVar1 - (var1Slot0*m_var1StepSize + m_var1Min))
     / m_var1StepSize);

  // Linear interpolation of the normalized value
  result
    = m_LUT[var1Slot0]*(1-var1Normalised)
    + m_LUT[var1Slot0+1]*var1Normalised;
}

void PCacheIAPWSIF97OneVar::init()
{
  m_properties.setClassName("PCacheIAPWSIF97OneVar");
  m_properties.setName("PCacheIAPWSIF97OneVar");
  m_properties.setDescription
    ("Computes a thermodynamic output variable or a transport "
     "coefficient 'out' from a thermodynamic variable 'var1', and "
     "possibly additional variables like 'var2' or constraints "
     "representing (partially) a thermodynamic state. For the "
     "specific meaning of 'out' and 'var1', and for the availability "
     "of additional thermodynamic variables or constraints for this "
     "specific module, see the "
     "end of this description (before the description of the "
     "attributes).\n"
     "Calculations are performed based on IAPWS-IF97 (International "
     "Association for the Properties of Water and Steam. Revised "
     "release on the IAPWS industrial formulation 1997 for the "
     "thermodynamic properties of water and steam. adadad, August "
     "2007) using the freesteam-library "
     "(http://freesteam.sourceforge.net/). Calculations are performed "
     "based on a preconstructed lookup table (LUT) in a predefined "
     "range, e.g. in 2D:\n"
     "[m_var1Min, m_var2Min] x [m_var1Max, m_var2Max] that is "
     "constructed once in advance to speedup the computation."
     ); 

  STRINGPC
      (symbol, m_symbolName,
       "Symbol name for the computed output 'out'.");
  STRINGPC
      (var1, m_var1Name,
       "Symbol name for 'var1'.");
  DOUBLEPC
      (var1Max, m_var1Max, 0,
       "Upper limit for the admissible range of 'var1'.");
  DOUBLEPC
      (var1Min, m_var1Min, 0,
       "Lower limit for the admissible range of 'var1'.");
  INTPC
      (LUTsizeVar1, m_arraySizeVar1, 0,
       "The number of support values for 'var1' in the generated "
       "look-up table.");
  
  m_var1Name = "undefined";
  m_var1Min = HUGE_VAL;
  m_var1Max = -HUGE_VAL;
  m_arraySizeVar1 = 0;
}

void PCacheIAPWSIF97OneVar::setup()
{

  if(m_var1Name == "undefined")
    throw gError("PCacheIAPWSIF97OneVar::setup for module " + className(), "Attribute 'var1' has value \"undefined\""); 
  if(m_var1Min == HUGE_VAL)
    throw gError("PCacheIAPWSIF97OneVar::setup for module " + className(), "Attribute 'var1Min' was not defined");
  if(m_var1Max == -HUGE_VAL)
    throw gError("PCacheIAPWSIF97OneVar::setup for module " + className(), "Attribute 'var1Max' was not defined");
  
  if(m_arraySizeVar1 == 0)
    throw gError("PCacheIAPWSIF97OneVar::setup for module " + className(), "Attribute 'LUTsizeVar1' was not defined.");

  if(m_var1Min >= m_var1Max)
    throw gError("PCacheIAPWSIF97OneVar::setup for module " + className(), "Attribute 'var1Min' (= " + ObjToString(m_var1Min) + ") >= attribute 'var1Max' (= " + ObjToString(m_var1Max) + ") not meaningful!");

  ParticleCache::setup();

  // setup the look-up table
  setupLUT();  

  // Now, if m_species == "all", we can freely copy this PCA including
  // m_LUT in function copyMySelf() called by registerCache()
  registerCache();

}


void PCacheIAPWSIF97OneVar::checkInputSymbolExistences(size_t colour) {
  
  if(Particle::s_tag_format[colour].attrExists(m_var1Name)) {
    if(m_datatype != Particle::s_tag_format[colour].attrByName(m_var1Name).datatype)
      throw gError("PressureCalculation::setup", "Var1 " + m_var1Name + " already exists as a non-scalar.");
    else m_var1Offset = Particle::s_tag_format[colour].offsetByName(m_var1Name);
  } else 
    throw gError("PressureCalculation::setup", "Symbol '" + m_var1Name + "' does not exist but required by this module.");

}


void PCacheIAPWSIF97OneVar::registerWithParticle()
{
}

