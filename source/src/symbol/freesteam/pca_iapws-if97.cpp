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


#ifdef HAVE_FREESTEAM

#include "pca_iapws-if97.h"
#include "simulation.h"
#include "manager_cell.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()

PCacheIAPWSIF97::PCacheIAPWSIF97
  (size_t colour, size_t offset, string symbolName)
    : ParticleCache(colour, offset, symbolName), m_LUT(NULL) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
}

PCacheIAPWSIF97::PCacheIAPWSIF97
    (/*Node*/Simulation* parent)
      : ParticleCache(parent), m_LUT(NULL) {
  m_stage = 0;
  m_datatype = DataFormat::DOUBLE;
  init();
}


PCacheIAPWSIF97::~PCacheIAPWSIF97()
{
  if(m_LUT) {
    delete [] m_LUT;
  }
}

void PCacheIAPWSIF97::setupLUT() {

  // auxiliary variables
  double deltaVar1 = 0.;
  double deltaVar2 = 0.;
  size_t slot;
  
  checkConstraints();
  
  // Step sizes depending on the size of the Array and the ranges of the input values.
  m_var1StepSize = (m_var1Max - m_var1Min)/(m_arraySizeVar1-1);
  m_var2StepSize = (m_var2Max - m_var2Min)/(m_arraySizeVar2-1);

  // Allocation of the 2D Array, depending on the given sizes.
  m_LUT = new double[m_arraySizeVar1*m_arraySizeVar2];

  // Calculate results at support points and store them into the LUT
  for(size_t i = 0; i < m_arraySizeVar1; ++i) {
    slot = i*m_arraySizeVar2;
    for (size_t j = 0; j < m_arraySizeVar2; ++j) {

      freesteamCalculationForState(m_LUT[slot],
				   m_var1Min + deltaVar1,
				   m_var2Min + deltaVar2);
      
      deltaVar2 += m_var2StepSize;
      ++slot;
    }
    deltaVar2 = 0.;
    deltaVar1 += m_var1StepSize; 
  }
}

void PCacheIAPWSIF97::calculateResult(double& result, const double& inputVar1, const double& inputVar2) const {

  // out of bounds?
  if((inputVar1 < m_var1Min) || (inputVar2 < m_var2Min)
     || (inputVar1 > m_var1Max) || (inputVar2 > m_var2Max))    
    throw gError("PCacheIAPWSIF97::calculateResult for module " + className(), "(" + m_var1Name + "," + m_var2Name + ") pair out of bounds: " + m_var1Name + "=" + ObjToString(inputVar1) + ", " + m_var2Name + "=" + ObjToString(inputVar2) + ", admissible [" + m_var1Name + "Min," + m_var1Name + "Max] = [" + ObjToString(m_var1Min) + "," + ObjToString(m_var1Max) + "], admissible [" + m_var2Name + "Min," + m_var2Name + "Max] = [" + ObjToString(m_var2Min) + "," + ObjToString(m_var2Max) + "].");

  // Calculation of the surrounding sampling points in 2D LUT
  size_t var1Slot0 = (floor((inputVar1-m_var1Min)/m_var1StepSize));
  size_t var2Slot0 = (floor((inputVar2-m_var2Min)/m_var2StepSize));
  // size_t var1Slot1 = var1Slot0+1;
  // size_t var2Slot1 = var2Slot0+1;

  // Normalisation of the input values.
  double var1Normalised =
    ((inputVar1 - (var1Slot0*m_var1StepSize + m_var1Min))
     / m_var1StepSize);
  double var2Normalised =
    ((inputVar2 - (var2Slot0*m_var2StepSize + m_var2Min))
     / m_var2StepSize);

  // Bilinear interpolation of the normalized values.
  // FIXME: Cache optimisation would require, e.g., that the 4 entries of the used square occupy adjacent array entries. Currently, at least two times 2 values are adjacent, which might already be good enough. Also a 2x2 square-wise ordering will fail (cache-wise) in 50%(?) of the cases because any square, shifted by 1 index in any direction is a valid square. Maybe bigger squares? Investigate if you have too much time ;)
  size_t var1Slot0Shift = var1Slot0*m_arraySizeVar2;
  size_t var1Slot1Shift = var1Slot0Shift + m_arraySizeVar2;
  result
    = m_LUT[var2Slot0 + var1Slot0Shift]
    *(1-var1Normalised)*(1-var2Normalised)
    + m_LUT[var2Slot0 + var1Slot1Shift]
    *var1Normalised*(1-var2Normalised)
    + m_LUT[var2Slot0 + 1 + var1Slot0Shift]
    *(1-var1Normalised)*var2Normalised
    + m_LUT[var2Slot0 + 1 + var1Slot1Shift]
    *var2Normalised*var1Normalised;
}

void PCacheIAPWSIF97::init()
{
  m_properties.setClassName("PCacheIAPWSIF97");
  m_properties.setName("PCacheIAPWSIF97");
  m_properties.setDescription
    ("Computes a thermodynamic output variable or a transport "
     "coefficient 'out' from two thermodynamic input variables 'var1' "
     "and 'var2' representing a thermodynamic state (var1, var2). For "
     "the specific "
     "meaning of 'out', 'var1', 'var2', see the end of this "
     "description (before the description of the attributes). "
     "Calculations are performed based on IAPWS-IF97 (International "
     "Association for the Properties of Water and Steam. Revised "
     "release on the IAPWS industrial formulation 1997 for the "
     "thermodynamic properties of water and steam. adadad, August "
     "2007) using the freesteam-library "
     "(http://freesteam.sourceforge.net/). Calculations are performed "
     "based on a preconstructed lookup table (LUT) in a predefined "
     "range [m_var1Min, m_var2Min] x [m_var1Max, m_var2Max] that is "
     "constructed once in advance to speedup the computation."
     ); 

  STRINGPC
      (symbol, m_symbolName,
       "Symbol name for the computed output 'out'.");
  STRINGPC
      (var1, m_var1Name,
       "Symbol name for 'var1'.");
  STRINGPC
      (var2, m_var2Name,
       "Symbol name for 'var2'.");
  DOUBLEPC
      (var1Max, m_var1Max, 0,
       "Upper limit for the admissible range of 'var1'.");
  DOUBLEPC
      (var1Min, m_var1Min, 0,
       "Lower limit for the admissible range of 'var1'.");
  DOUBLEPC
      (var2Max, m_var2Max, 0,
       "Upper limit for the admissible range of 'var2'.");
  DOUBLEPC
      (var2Min, m_var2Min, 0,
       "Lower limit for the admissible range of 'var2'.");
  INTPC
      (LUTsizeVar1, m_arraySizeVar1, 0,
       "The number of support values for 'var1' in the generated "
       "look-up table.");
  INTPC
      (LUTsizeVar2, m_arraySizeVar2, 0,
       "The number of support values for 'var2' in the generated "
       "look-up table.");
  
  m_var1Name = "undefined";
  m_var2Name = "undefined";
  m_var1Min = HUGE_VAL;
  m_var2Min = HUGE_VAL;
  m_var1Max = -HUGE_VAL;
  m_var2Max = -HUGE_VAL;
  m_arraySizeVar1 = 0;
  m_arraySizeVar2 = 0;
}

void PCacheIAPWSIF97::setup()
{

  if(m_var1Name == "undefined")
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var1' has value \"undefined\""); 
  if(m_var2Name == "undefined")
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var2' has value \"undefined\""); 

  if(m_var1Min == HUGE_VAL)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var1Min' was not defined");
  if(m_var2Min == HUGE_VAL)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var2Min' was not defined");
  if(m_var1Max == -HUGE_VAL)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var1Max' was not defined");
  if(m_var2Max == -HUGE_VAL)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var2Max' was not defined");
  
  if(m_arraySizeVar1 == 0)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'LUTsizeVar1' was not defined.");
  if(m_arraySizeVar2 == 0)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'LUTsizeVar2' was not defined.");

  if(m_var1Min >= m_var1Max)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var1Min' (= " + ObjToString(m_var1Min) + ") >= attribute 'var1Max' (= " + ObjToString(m_var1Max) + ") not meaningful!");
  if(m_var2Min >= m_var2Max)
    throw gError("PCacheIAPWSIF97::setup for module " + className(), "Attribute 'var2Min' (= " + ObjToString(m_var2Min) + ") >= attribute 'var1Max' (= " + ObjToString(m_var2Max) + ") not meaningful!");

  ParticleCache::setup();
  
  // setup the look-up table
  setupLUT();  

  // Now, if m_species == "all", we can freely copy this PCA including
  // m_LUT in function copyMySelf() called by registerCache()
  registerCache();

}


void PCacheIAPWSIF97::checkInputSymbolExistences(size_t colour) {

  if(Particle::s_tag_format[colour].attrExists(m_var2Name)) {
    if(m_datatype != Particle::s_tag_format[colour].attrByName(m_var2Name).datatype)
      throw gError("PressureCalculation::setup", "Var2 " + m_var2Name + " already exists as a non-scalar.");
    else m_var2Offset = Particle::s_tag_format[colour].offsetByName(m_var2Name);
  } else 
    throw gError("PressureCalculation::setup", "Symbol '" + m_var2Name + "' does not exist but required by this module.");
  
  if(Particle::s_tag_format[colour].attrExists(m_var1Name)) {
    if(m_datatype != Particle::s_tag_format[colour].attrByName(m_var1Name).datatype)
      throw gError("PressureCalculation::setup", "Var1 " + m_var1Name + " already exists as a non-scalar.");
    else m_var1Offset = Particle::s_tag_format[colour].offsetByName(m_var1Name);
  } else 
    throw gError("PressureCalculation::setup", "Symbol '" + m_var1Name + "' does not exist but required by this module.");

}


void PCacheIAPWSIF97::registerWithParticle()
{
}

#endif    // HAVE_FREESTEAM

