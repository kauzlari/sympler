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


#include "config.h"

#include "vtk_support.h"


/*--- ASCII writers ---*/

void vtkWriteVectorInt_Ascii(ostream &s, string name, vector_int_sp vi)
{
    //    MSG_DEBUG("vtkWriteVectorDouble", name);

    s << "SCALARS " << name << " int 1" << endl;
    s << "LOOKUP_TABLE default" << endl;

    for (vector<int>::iterator i = vi->begin(); i != vi->end(); i++) {
        s << *i << endl;
    }

    s << endl;
}

void vtkWriteVectorDouble_Ascii(ostream &s, string name, vector_double_sp vd)
{
    //    MSG_DEBUG("vtkWriteVectorDouble", name);

    s << "SCALARS " << name << " double 1" << endl;
    s << "LOOKUP_TABLE default" << endl;

    for (vector<double>::iterator i = vd->begin(); i != vd->end(); i++) {
        s << *i << endl;
    }

    s << endl;
}


void vtkWriteVectorPoint_Ascii(ostream &s, string name, vector_point_sp vp)
{
    //    MSG_DEBUG("vtkWriteVectorPoint", name);

    s << "VECTORS " << name << " double" << endl;

    for (vector<point_t>::iterator i = vp->begin(); i != vp->end(); i++) {
        s << i->x << " " << i->y << " " << i->z << endl;
    }

    s << endl;
}


void vtkWriteVectorTensor_Ascii(ostream &s, string name, vector_tensor_sp vp)
{
    s << "TENSORS " << name << " double" << endl;

    for (vector<tensor_t>::iterator i = vp->begin(); i != vp->end(); i++) {
      for (int a = 0; a < SPACE_DIMS; ++a) {
        for (int b = 0; b < SPACE_DIMS; ++b)
          s << (*i)(a, b) << " ";
        s << endl;
      }
      s << endl;
    }

    s << endl;
}


/*--- Binary writers ---*/

void vtkWriteVectorInt_Binary(ostream &s, string name, vector_int_sp vi)
{
    //    MSG_DEBUG("vtkWriteVectorDouble", name);

    s << "SCALARS " << name << " int 1" << endl;
    s << "LOOKUP_TABLE default" << endl;

    for (vector<int>::iterator i = vi->begin(); i != vi->end(); i++) {
        int num = (*i);
        vtkSwapInt(num);
        s.write((char*) &num, sizeof(int));
    }

    s << endl;
}

void vtkWriteVectorDouble_Binary(ostream &s, string name, vector_double_sp vd)
{
    //    MSG_DEBUG("vtkWriteVectorDouble", name);

    s << "SCALARS " << name << " double 1" << endl;
    s << "LOOKUP_TABLE default" << endl;

    for (vector<double>::iterator i = vd->begin(); i != vd->end(); i++) {
        double d = *i;
        vtkSwapDouble(d);
        s.write((char*) &d, sizeof(double));
    }

    s << endl;
}


void vtkWriteVectorPoint_Binary(ostream &s, string name, vector_point_sp vp)
{
    //    MSG_DEBUG("vtkWriteVectorPoint", name);

    s << "VECTORS " << name << " double" << endl;

    for (vector<point_t>::iterator i = vp->begin(); i != vp->end(); i++) {
        point_t p = *i;
        vtkSwapPoint(p);
        s.write((char*) &p.x, sizeof(double));
        s.write((char*) &p.y, sizeof(double));
        s.write((char*) &p.z, sizeof(double));
    }

    s << endl;
}


void vtkWriteVectorTensor_Binary(ostream &s, string name, vector_tensor_sp vp)
{
    s << "TENSORS " << name << " double" << endl;

    for (vector<tensor_t>::iterator i = vp->begin(); i != vp->end(); i++) {
      for (int a = 0; a < SPACE_DIMS; ++a) {
        for (int b = 0; b < SPACE_DIMS; ++b)
        {
          double d = (*i)(a, b);
          vtkSwapDouble(d);
          s.write((char*) &d, sizeof(double));
          //           s.write((char*) &(*i)(a, b), sizeof(double));
      
        }
      }
    }

    s << endl;
}


/*--- Byte swappers ---*/


#ifdef IS_BIG_ENDIAN

void vtkSwapDouble(double &data) { }
void vtkSwapInt(int &data) { }
void vtkSwapPoint(point_t &data) { }

#else

typedef void (*swap_t)(char*);


void vtkSwap2Bytes(char* data)
{
    char one_byte;  
    one_byte = data[0]; data[0] = data[1]; data[1] = one_byte;
}


void vtkSwap4Bytes(char* data)
{ 
    char one_byte; 
    one_byte = data[0]; data[0] = data[3]; data[3] = one_byte;
    one_byte = data[1]; data[1] = data[2]; data[2] = one_byte; 
}


void vtkSwap8Bytes(char* data)
{ 
  char one_byte;
  one_byte = data[0]; data[0] = data[7]; data[7] = one_byte;
  one_byte = data[1]; data[1] = data[6]; data[6] = one_byte;
  one_byte = data[2]; data[2] = data[5]; data[5] = one_byte;
  one_byte = data[3]; data[3] = data[4]; data[4] = one_byte; 
}


swap_t __swappers[8] = {
    NULL, vtkSwap2Bytes, NULL, vtkSwap4Bytes, 
    NULL, NULL, NULL, vtkSwap8Bytes
};


void vtkSwapDouble(double &data)
{
    assert(__swappers[sizeof(double)-1] != NULL);
    __swappers[sizeof(double)-1]((char*) &data);
}


void vtkSwapInt(int &data)
{
    assert(__swappers[sizeof(int)-1] != NULL);
    __swappers[sizeof(int)-1]((char*) &data);
}


void vtkSwapPoint(point_t &data)
{
    for (int i = 0; i < SPACE_DIMS; i++)
        vtkSwapDouble(data[i]);
}


#endif

