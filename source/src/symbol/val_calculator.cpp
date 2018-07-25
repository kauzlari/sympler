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



#include "manager_cell.h"
#include "val_calculator.h"
#include "simulation.h"
#include "colour_pair.h"

#define M_SIMULATION ((Simulation *) m_parent)
#define M_PHASE M_SIMULATION->phase()
#define M_MANAGER M_PHASE->manager()
#define M_CONTROLLER M_SIMULATION->controller()

using namespace std;

ValCalculator::ValCalculator(string symbol): Symbol(symbol)
{
  m_datatype = DataFormat::DOUBLE;
}

ValCalculator::ValCalculator(/*Node*/Simulation* parent): Symbol(parent)
{
  m_datatype = DataFormat::DOUBLE;
  init();
}

void ValCalculator::init()
{
  m_properties.setClassName("ValCalculator");

  STRINGPC
      (species1, m_species.first,
       "Name for the first species of the pairs, this Symbol is used for.");
  
  STRINGPC
      (species2, m_species.second,
       "Name for the second species of the pairs, this Symbol is used for.");
  
  m_species.first = "undefined";
  m_species.second = "undefined";
 
}

void ValCalculator::cleanSymbol(string& name) const
{
  if(name[0] == '{' || name[0] == '[') {
    // remove the first bracket
    name.erase(0, 1);
    // remove the last bracket; don't know why, but with these arguments it works
    name.erase(name.size()-1, name.size()-1);
  }
  // and the "i", "j" and "ij" of the pair expression have to be removed
  if(name[name.size()-2] == 'i' && name[name.size()-1] == 'j')
    // remove "ij"
    name.erase(name.size()-2, name.size()-1); 
  else
    // remove "i" or "j"
    name.erase(name.size()-1, name.size()-1);
  MSG_DEBUG("ValCalculator::cleanPairSymbol", className() << ": shortened name of symbol: " << name);
}



