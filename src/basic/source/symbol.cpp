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


#include "symbol.h"
#include "node.h"

using namespace std;

REGISTER_SMART_ENUM
(SymbolFactory,
 "Expression, computed and saved in a symbol, to be used in expression of other modules.");

Symbol::Symbol(string symbol)
  : Node(NULL), m_stage(-1), m_phase(1), m_symbolName(symbol/*"undefined"*/), m_datatype(DataFormat::DOUBLE), m_overwrite(false)
{
}

Symbol::Symbol(/*Node*/Simulation *parent)
  : Node((Node*) parent), m_stage(-1), m_phase(1)
{
  init();
}
 
/* Methods */

void Symbol::init()
{
  m_properties.setClassName("Symbol");

  INTPC
      (stage, m_phaseUser, -1,
       "When during the timestep should this Calculator be called? Current possibilities are: \n\"0\": At the very beginning of the timestep.\n\"1\": Right before the force calculations.\n\"2\": At both instances.");

//   BOOLPC
//       (overwrite, m_overwrite,
//        "Is this calculator allowed to overwrite already existing symbols " 
//            "with name 'symbol' ?");

  m_symbolName = "undefined";
  m_phaseUser = 1;
  // We set this as the default value for m_overwrite even though only selected 
  // sub-classes provide an XML-attribute for this member
  m_overwrite = false;
}


#ifdef _OPENMP
int Symbol::setNumOfDoubles() {
    return DataFormat::getNumOfDoubles(m_datatype);
}
#endif
