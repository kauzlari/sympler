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


#include "function_arbitrary.h"

#include "fp_scalar.h"
#include "fp_vector.h"
#include "fp_tensor.h"

#include "weighting_function.h"

/* Function */

FunctionArbitrary::FunctionArbitrary()
  : m_returnType(Variant::SCALAR)/*, m_default_wf(NULL)*/
{
}


FunctionArbitrary::~FunctionArbitrary()
{
}


string FunctionArbitrary::int2CExpression(string base_name, double base_offset)
{
  // the outermost bracket is needed, e.g., to avoid the merging of "/" and 
  // "*" to "/*" !!!
  return
      "(*((int*) ((char*) " + base_name + " + " +
      ObjToString(base_offset) +
      ")))";
}

string FunctionArbitrary::double2CExpression(string base_name, double base_offset)
{
  // the outermost bracket is needed, e.g., to avoid the merging of "/" and 
  // "*" to "/*" !!!
  return
      "(*((double*) ((char*) " + base_name + " + " +
      ObjToString(base_offset) +
      ")))";
}


void FunctionArbitrary::point2CExpression(string *str_vec, string base_name, double base_offset)
{
  // the outermost bracket is needed, e.g., to avoid the merging of "/" and 
  // "*" to "/*" !!!
  for (int i = 0; i < SPACE_DIMS; ++i) {
    size_t offset = 0;

    switch (i) {
      case 0: offset = offsetof(point_t, x); break;
      case 1: offset = offsetof(point_t, y); break;
      case 2: offset = offsetof(point_t, z); break;
    }

    str_vec[i] =
        "(*((double*) ((char*) " + base_name + " + " +
        ObjToString(base_offset + offset) +
        ")))";
  }
}

void FunctionArbitrary::tensor2CExpression(string str_vec[SPACE_DIMS][SPACE_DIMS], string base_name, double base_offset)
{
  for (int i = 0; i < 9; ++i) {
    size_t offset = 0;

    switch (i) {
      case 0: 
      {
        offset = offsetof(tensor_t, xx);
        str_vec[0][0] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 1: 
      {
        offset = offsetof(tensor_t, xy);
        str_vec[0][1] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 2: 
      {
        offset = offsetof(tensor_t, xz);
        str_vec[0][2] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 3: 
      {
        offset = offsetof(tensor_t, yx);
        str_vec[1][0] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 4: 
      {
        offset = offsetof(tensor_t, yy);
        str_vec[1][1] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 5: 
      {
        offset = offsetof(tensor_t, yz);
        str_vec[1][2] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
//         MSG_DEBUG("FunctionArbitrary::tensor2CExpression", "yz-case: offset = " << offset << "string = " << str_vec[1][2]);
      } break;
        
      case 6: 
      {
        offset = offsetof(tensor_t, zx);
        str_vec[2][0] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 7: 
      {
        offset = offsetof(tensor_t, zy);
        str_vec[2][1] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
      case 8: 
      {
        offset = offsetof(tensor_t, zz);
        str_vec[2][2] =
            "(*((double*) ((char*) " + base_name + " + " +
            ObjToString(base_offset + offset) +
            ")))";
      } break;
        
/*      case 1: offset = offsetof(tensor_t, xy); break;
      case 2: offset = offsetof(tensor_t, xz); break;
      case 3: offset = offsetof(tensor_t, yx); break;
      case 4: offset = offsetof(tensor_t, yy); break;
      case 5: offset = offsetof(tensor_t, yz); break;
      case 6: offset = offsetof(tensor_t, zx); break;
      case 7: offset = offsetof(tensor_t, zy); break;
      case 8: offset = offsetof(tensor_t, zz); break;*/
    }

  // the outermost bracket is needed, e.g., to avoid the merging of "/" and 
  // "*" to "/*" !!!
/*    str_vec[i] =
    "(*((double*) ((char*) " + base_name + " + " +
    ObjToString(base_offset + offset) +
    ")))";*/
  }
}


void FunctionArbitrary::addInt(string name, string base_name, size_t offset)
{
  m_parser.addSymbol
      (new FPScalarVariable
      (name,
       NULL,
       int2CExpression(base_name, offset)));
}


void FunctionArbitrary::addDouble(string name, string base_name, size_t offset)
{
  m_parser.addSymbol
      (new FPScalarVariable
      (name,
       NULL,
       double2CExpression(base_name, offset)));
}


void FunctionArbitrary::addPoint(string name, string base_name, size_t offset)
{
  string str_vec[SPACE_DIMS];

  point2CExpression(str_vec, base_name, offset);
  m_parser.addSymbol
      (new FPVectorVariable
      ("[" + name + "]",
       NULL,
       str_vec));
}


void FunctionArbitrary::addTensor(string name, string base_name, size_t offset)
{
  string str_tensor[SPACE_DIMS][SPACE_DIMS];

  tensor2CExpression(str_tensor, base_name, offset);
  m_parser.addSymbol
      (new FPTensorVariable
      ("{" + name + "}",
       NULL,
       str_tensor));
}


string FunctionArbitrary::addSuffix(string symbol, string suffix)
{
  string s(symbol);
  int pos;

  pos = s.find('[');

  if (pos != -1)
    s.insert(pos, suffix);
  else
    s.append(suffix);
  
  return s;
}


void FunctionArbitrary::addAllFromDataFormat(string base_name, DataFormat *format, string suffix)
{
  for (size_t i = 0; i < format->rows(); ++i) {
    DataFormat::attribute_t attr = format->attrByIndex(i);

    switch (attr.datatype) {
      case DataFormat::DOUBLE:
      {
//       MSG_DEBUG("Function::addAllFromDataFormat", "adding symbol " << attr.symbol << " with suffix " << suffix);
        addDouble(addSuffix(attr.symbol, suffix), base_name, attr.offset);
      }
      break;
      case DataFormat::POINT:
        addPoint(addSuffix(attr.symbol, suffix), base_name, attr.offset);
        break;
      case DataFormat::TENSOR:
        addTensor(addSuffix(attr.symbol, suffix), base_name, attr.offset);
        break;
      default:
      // for everything else, nothing is added
      // the throwing of an exception was commented out (2006/01/10)
//       throw gError
        MSG_DEBUG
            ("Function::addAllFromDataFormat",
             "Don't know how to handle data type " + ObjToString(attr.datatype) +
                 " for entry '" + attr.name + "'");
    }
  }
}

void FunctionArbitrary::compile()
{
  vector<string> retstr;

  switch (m_returnType) {
    case Variant::SCALAR: {
      retstr.push_back(double2CExpression("result", 0));
    }
    break;
    case Variant::VECTOR: {
      point_t p;

      for (int i = 0; i < SPACE_DIMS; ++i) {
      //      retstr.push_back(double2CExpression("result", offsetof(point_t, coords[i])));
        retstr.push_back(double2CExpression("result", (size_t) &p.coords[i] - (size_t) &p));
      }
    }
    break;
    case Variant::TENSOR: {
      tensor_t t;

      for (int i = 0; i < SPACE_DIMS; ++i) 
        for (int j = 0; j < SPACE_DIMS; ++j) {
// 		retstr.push_back(double2CExpression("result", offsetof(tensor_t, tensor[i][j])));
        retstr.push_back(double2CExpression("result", (size_t) &t.tensor[j+i*SPACE_DIMS]/*[i][j]*/ - (size_t) &t));
        }
    }
    break;
    default:
      throw gError
          ("FunctionParticle::compile",
           "Don't know how to handle return type " + ObjToString(m_returnType));
  }

  m_compiler.setResultStrings(retstr);
}

// void FunctionArbitrary::setWeightF(WeightingFunction *wf) {
//   m_default_wf = wf;
// }

// Function::Function(const Function&)
// {
//   MSG_DEBUG("Function::Function(const Function&)", "COPY CONSTRUCTOR called"); 
// }
