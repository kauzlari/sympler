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


#include "unary_operators.h"

#include "typed_value.h"

/* UnaryFunction */

FNUnaryOperator::FNUnaryOperator(string name): FunctionNode(name), m_a(NULL)
{
}


FNUnaryOperator::~FNUnaryOperator()
{
  if (!dynamic_cast<TypedValue*>(m_a))
  {
//     MSG_DEBUG("FNUnaryOperator::~FNUnaryOperator", "really DELETING");

    delete m_a;

  }
}

void FNUnaryOperator::setUnary(FunctionNode *a)
{
  if (!dynamic_cast<TypedValue*>(m_a))
    delete m_a;

  m_a = a;
}



FNNegation::FNNegation(): FNUnaryOperator("-")
{
}							


FNNegation::~FNNegation()
{				
}							

      							
Variant FNNegation::value() const
{	
  Variant::variant_type_t t;					
  Variant vca = m_a->value();					
									
  t = vca.typeId();							
									
  if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) { 
    Variant r(t);							
    vector_double_t &vr = r.doubles();				
    vector_double_t va = vca.doubles();				
      
    for (size_t i = 0; i < va.size(); ++i) {				
      vr[i] = -va[i];
    }								
      
    return r;							
  } else								
    throw gError							
      ("FNNegation::value",						
       "Can't handle Variant type id " + ObjToString(t) + ".");	
}									
    									
Variant FNNegation::toC() const
{					
  Variant::variant_type_t t;					
  Variant vca = m_a->toC();						
									
  t = vca.typeId();							
    
  if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	
      t == Variant::TENSOR_STRING) {				
    Variant r(t);							
    vector_string_t &vr = r.strings();				
    vector_string_t va = vca.strings();				
									
    for (size_t i = 0; i < va.size(); ++i) {				
      vr[i] = "(-" + va[i] + ")";
    }								
									
    return r;							
  } else								
    throw gError							
      ("FNNegation::toC",						
       "Can't handle Variant type id " + ObjToString(t) + ".");	
}
									
//const FunctionNode_Register<FNNegation> fn_step(FunctionNode::UNARY_OPERATOR, 6, "-");
