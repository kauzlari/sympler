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


#ifndef __FUNCTION_ARBITRARY_H
#define __FUNCTION_ARBITRARY_H 

using namespace std;
/* needed for offsetof(type, member) */
#include <cstddef>
#include "function.h"
#include "variant.h"
#include "data_format.h"


class WeightingFunction;

/*!
 * Base class for function which takes structures and/or data
 * as their arguments
 */
class FunctionArbitrary : public Function
{
  protected:
  /*!
     * Return type
   */
    Variant::variant_type_t m_returnType;



  /*!
     * Add an integer parameter of name \a name. \a base_name defines the structure
     * relative of which (at offset \a offset) to find the double.
     * @param name The name of this double (in the expression)
     * @param base_name Pointer to the structure that contains this double (in the C file)
     * @param offset Offset relative to the strucutre where to find this double (in the C file)
   */
    void addInt(string name, string base_name, size_t offset);

  /*!
     * Add a double parameter of name \a name. \a base_name defines the structure
     * relative of which (at offset \a offset) to find the double.
     * @param name The name of this double (in the expression)
     * @param base_name Pointer to the structure that contains this double (in the C file)
     * @param offset Offset relative to the strucutre where to find this double (in the C file)
   */
    void addDouble(string name, string base_name, size_t offset);

  /*!
     * Add a point parameter of name \a name. \a base_name defines the structure
     * relative from which (at offset \a offset) to find the point.
     * @param name The name of this point (in the expression)
     * @param base_name Pointer to the structure that contains this point (in the C file)
     * @param offset Offset relative to the structure where to find this point (in the C file)
   */
    void addPoint(string name, string base_name, size_t offset);

  /*!
     * Add a tensor parameter of name \a name. \a base_name defines the structure
     * relative from which (at offset \a offset) to find the tensor.
     * @param name The name of this tensor (in the expression)
     * @param base_name Pointer to the structure that contains this tensor (in the C file)
     * @param offset Offset relative to the structure where to find this tensor (in the C file)
   */
    void addTensor(string name, string base_name, size_t offset);

  /*!
     * Add all attributes from \a DataFormat \a format. The \a suffix is appended to
     * the short description given in the \a DataFormat.
     * @param base_name Pointer to the tag that contains this information (in the C file)
     * @param format Format description for the tag
     * @param suffix Suffix to append to all attribute names found in \a format
   */
    void addAllFromDataFormat(string base_name, DataFormat *format, string suffix = "");
                                        
  public:
  /*!
   * Constructor
   */
    FunctionArbitrary();

  /*!
     * Destructor
   */
    virtual ~FunctionArbitrary();
  
  /*!
     * Compile this function using gcc
   */
    virtual void compile();

  /*!
     * Set the return type of this function, i.e., scalar, vector or tensor
   */
    void setReturnType(Variant::variant_type_t type) {
      m_returnType = type;
    }

    Variant::variant_type_t returnType()
    {
      return m_returnType; 
    }
    
  
  /*!
     * Return the appropriate expression in the C file for access to an int argument
     * @param base_name Pointer to the base structure
     * @param base_offset Offset relative to the base structure
   */
    static string int2CExpression(string base_name, double base_offset);
  
  /*!
     * Return the appropriate expression in the C file for access to a double argument
     * @param base_name Pointer to the base structure
     * @param base_offset Offset relative to the base structure
   */
    static string double2CExpression(string base_name, double base_offset);

  /*!
     * Return the appropriate expression in the C file for access to a point argument
     * @param str_vec String for each vector entry (return value)
     * @param base_name Pointer to the base structure
     * @param base_offset Offset relative to the base structure
   */
    static void point2CExpression(string *str_vec, string base_name, double base_offset);

  /*!
     * Return the appropriate expression in the C file for access to a tensor argument
     * @param str_vec String for each tensor entry (return value)
     * @param base_name Pointer to the base structure
     * @param base_offset Offset relative to the base structure
   */
    static void tensor2CExpression(string str_vec[SPACE_DIMS][SPACE_DIMS], string base_name, double base_offset);

  /*!
     * Add a suffix to symbol \a symbol, which mean if there is a [
     * present add the suffix in front of the [
     * @param symbol Symbol to add the suffix to
     * @param suffix Suffix to add
   */
    static string addSuffix(string symbol, string suffix);
};

#endif
