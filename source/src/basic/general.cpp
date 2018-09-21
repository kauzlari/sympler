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



#include "general.h"

#ifdef _OPENMP
  #include "omp.h"
#endif

#include <iomanip>

using namespace std;

double global::R = 8.31441;

/*!
 * make_filename takes the filename 's' and appends the number 'n'
 */
string make_filename(string s, int n)
{
  // if there is a dot, place the number n behind it (rfind searches from the end)
  int p = s.rfind('.');
  // otherwise put the number n at the end
  if(p<0)
    p = s.size();
  stringstream h;
  
  h << setfill('0') << setw(5) << n;
  
  s.insert(p, "_"+h.str());
  
  return s;
}


bool g_stringIsInPipeList(string s, string pipeStringList)
{

  if (pipeStringList != "---") {
    bool run = true;
    string working = pipeStringList;
    while(run) {
      string cur;
      size_t pos = working.find('|');

      if (pos == string::npos) {
      	run = false;
      	cur = working;
      }
      else {
      	cur = string(working, 0, pos);
      	working = string(working, pos+1);
      }

      if(s == cur)
      	return true; // should quit the while loop
    }

  }

  return false;
}


