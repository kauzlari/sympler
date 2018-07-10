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


#include <assert.h>

#include <string>
#include <fstream>
#include <iostream>

#include "vtk_support.h"

using namespace std;


void vtk2profile(string vtkfn, string data_name)
{
  static const string cell_data = "CELL_DATA";
  static const string token1_scalar = "SCALARS";
  static const string token1_vector = "VECTORS";
  static const string token2 = "LOOKUP_TABLE";

  ifstream in;
  int size, i, n_data;
  char *buffer;
  bool tok1_found, cd_found, data_found, data_is_scalar;

  in.open(vtkfn.c_str(), ifstream::binary);

  if (!in.good()) {
    cout << "Error opening '" << vtkfn << "'" << endl;
    exit(1);
  }

  in.seekg (0, ios::end);
  size = in.tellg();
  in.seekg (0, ios::beg);

  buffer = (char*) malloc(size);

  in.read(buffer, size);

  i = 0;
  cd_found = false;
  tok1_found = false;
  data_found = false;
  data_is_scalar = true;
  while (i < size) {
    if (buffer[i] == cell_data[0] && !cd_found && !tok1_found && !data_found) {
      if (!memcmp(&buffer[i], cell_data.c_str(), cell_data.size())) {
        int j;
        char s[80];

        cd_found = true;
        i += cell_data.size();

        while (buffer[i] == ' ')
          i++;

        j = 0;
        while (buffer[i+j] != '\n')
          j++;

        if (j >= 80) {
          cout << "j >= 80" << endl;
          exit(1);
        }

        memcpy(s, &buffer[i], j);
        s[j] = 0;

        n_data = atoi(s);
      } else
        i++;
    } else if ((buffer[i] == token1_scalar[0] || buffer[i] == token1_vector[0]) && cd_found && !tok1_found && !data_found) {
      if (!memcmp(&buffer[i], token1_scalar.c_str(), token1_scalar.size())) {
        /* Is this SCALARS? */
        i += token1_scalar.size();
        while (buffer[i] == ' ')
          i++;
        if (!memcmp(&buffer[i], data_name.c_str(), data_name.size())) {
          i += data_name.size();
          tok1_found = true;
          data_is_scalar = true;
        }
      } else if (!memcmp(&buffer[i], token1_vector.c_str(), token1_vector.size())) {
        /* Or VECTORS? */
        i += token1_vector.size();
        while (buffer[i] == ' ')
          i++;
        if (!memcmp(&buffer[i], data_name.c_str(), data_name.size())) {
          i += data_name.size();
          tok1_found = false;
          data_found = true;
          data_is_scalar = false;

          while (buffer[i] != '\n')
            i++;
          i++;

          for (int j = 0; j < n_data; j++) {
            double d1 = *((double*) (&buffer[i]));
            double d2 = *((double*) (&buffer[i+sizeof(double)]));
            double d3 = *((double*) (&buffer[i+2*sizeof(double)]));

            point_t p = { { { d1, d2, d3 } } };

            vtkSwapPoint(p);

            cout << p.x << " " << p.y << " " << p.z << endl;

            i += 3*sizeof(double);
          }
        }
      } else
        i++;
    } else if (buffer[i] == token2[0] && tok1_found && cd_found && !data_found && data_is_scalar) {
      if (!memcmp(&buffer[i], token2.c_str(), token2.size())) {
        data_found = true;
        /* Okay, here we are. */
        i += token2.size();
        while (buffer[i] != '\n')
          i++;

        i++;

        for (int j = 0; j < n_data; j++) {
          double d = *((double*) (&buffer[i]));
          
          vtkSwapDouble(d);

          cout << d << endl;

          i += sizeof(double);
        }
          
        tok1_found = false;
      } else
        i++;
    } else
      i++;
  }

  free(buffer);

  in.close();

  if (!data_found)
    cerr << "!!! cell data '" << data_name << "' not found !!!" << endl;
}



int main(int argc, char *argv[])
{
  if (argc != 3) {
    cout << "Syntax: vtk2profile <vtk-filename> <scalar-name>" << endl;
  } else {
    vtk2profile(argv[1], argv[2]);
  }
}
