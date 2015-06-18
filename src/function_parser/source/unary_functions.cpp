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


#include "unary_functions.h"

#include "typed_value.h"

/* UnaryFunction */

FNUnaryFunction::FNUnaryFunction(string name): FunctionNode(name), m_a(NULL)
{
}


FNUnaryFunction::~FNUnaryFunction()
{
  if (!dynamic_cast<TypedValue*>(m_a))
  {
//     MSG_DEBUG("FNUnaryFunction::~FNUnaryFunction", "really DELETING");
    delete m_a;
  }
}


void FNUnaryFunction::setUnary(FunctionNode *a)
{
  if (!dynamic_cast<TypedValue*>(m_a))
    delete m_a;

  m_a = a;
}


#define MAKE_UNARY_FUNCTION(str, expr, priority)	\
  /*! Definition of the str function */			\
  class FN ## str: public FNUnaryFunction		\
  {							\
  public:						\
    FN ## str(): FNUnaryFunction(#str) {		\
    }							\
      							\
    virtual ~FN ## str() {				\
    }							\
      							\
    virtual Variant value() const {					\
      Variant::variant_type_t t;					\
      Variant vca = m_a->value();					\
									\
      t = vca.typeId();							\
									\
      if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) { \
	Variant r(t);							\
	vector_double_t &vr = r.doubles();				\
	vector_double_t va = vca.doubles();				\
									\
	for (size_t i = 0; i < va.size(); ++i) {				\
	  vr[i] = expr(va[i]);						\
	}								\
									\
	return r;							\
      } else								\
	throw gError							\
	  ("FN" #str "::value",						\
	   "Can't handle Variant type id " + ObjToString(t) + ".");	\
    }									\
    									\
    virtual Variant toC() const {					\
      Variant::variant_type_t t;					\
      Variant vca = m_a->toC();						\
									\
      t = vca.typeId();							\
									\
      if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	\
	  t == Variant::TENSOR_STRING) {				\
	Variant r(t);							\
	vector_string_t &vr = r.strings();				\
	vector_string_t va = vca.strings();				\
									\
	for (size_t i = 0; i < va.size(); ++i) {				\
	  vr[i] = #expr "(" + va[i] + ")";				\
	}								\
									\
	return r;							\
      } else								\
	throw gError							\
	  ("FN" #str "::toC",						\
	   "Can't handle Variant type id " + ObjToString(t) + ".");	\
    }									\
  };									\
									\
  const FunctionNode_Register<FN ## str>				\
  fn_ ## str(FunctionNode::UNARY_FUNCTION, priority, #str)

MAKE_UNARY_FUNCTION(sqrt, sqrt, 6);

MAKE_UNARY_FUNCTION(sin, sin, 6);
MAKE_UNARY_FUNCTION(cos, cos, 6);
MAKE_UNARY_FUNCTION(tan, tan, 6);
MAKE_UNARY_FUNCTION(atan, atan, 6);

// next two only defined for arguments in [-1,1] !!!
// User must make sure this holds, otherwise the functions produce nan.
MAKE_UNARY_FUNCTION(asin, asin, 6);
MAKE_UNARY_FUNCTION(acos, acos, 6);

MAKE_UNARY_FUNCTION(sinh, sinh, 6);
MAKE_UNARY_FUNCTION(cosh, cosh, 6);
MAKE_UNARY_FUNCTION(tanh, tanh, 6);

MAKE_UNARY_FUNCTION(abs, fabs, 6);
MAKE_UNARY_FUNCTION(exp, exp, 6);
MAKE_UNARY_FUNCTION(round, round, 6);


class FNDet: public FNUnaryFunction	      
{							
  public:						
    FNDet(): FNUnaryFunction("det") {		
    }							

    virtual ~FNDet() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      // the return type is always a scalar, independently of the argument 
      Variant r(Variant::SCALAR);
      
      if (t != Variant::TENSOR) throw gError("FNDet::value", "The det() function may only be used with tensorial arguments");
      else
      {
        double& d = r.scalar();
        vector_double_t va = vca.doubles();				
        assert(va.size() == 9);
        
        d = va[2]*(va[3]*va[7] - va[4]*va[6]) 
	  + va[1]*(va[5]*va[6] - va[3]*va[8]) 
	  + va[0]*(va[4]*va[8] - va[5]*va[7]);
      }
      return r;
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::TENSOR_STRING) {				
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR_STRING);							
        // this one should have length one							
        vector_string_t va = vca.strings();				
        assert(va.size() == 9);
									
        r.scalarString() = 
	  "(" + va[2]+ "*" + "(" + va[3] + "*" + va[7] + "-" + va[4] + "*" + va[6] + ")" 
	   + "+" +  va[1] + "*" + "(" + va[5] + "*" + va[6] + "-" + va[3] + "*" + va[8] + ")" 
	   + "+" + va[0] + "*" + "(" + va[4] + "*" + va[8] + "-" + va[5] + "*" + va[7] + "))";

	return r;							
      } 
      else								
        throw gError							
            ("FNDet::toC",						
             "The det() function may only be used with tensorial arguments");	
    }									
};									
									
const FunctionNode_Register<FNDet> fn_det(FunctionNode::UNARY_FUNCTION, 6, "det");


class FNDiagT: public FNUnaryFunction	      
{							
  public:						
    FNDiagT(): FNUnaryFunction("diagMat") {		
    }							

    virtual ~FNDiagT() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::VECTOR) throw gError("FNDiagT::value", "The diagMat() function may only be used with vector arguments");
      else
      {
        // the return type is always a tensor 
        Variant vr(Variant::TENSOR);
        vector_double_t &r = vr.doubles();
        vector_double_t a = vca.doubles();
        // make sure that the off-diagonal elements are not undefined				
        for(size_t i = 0; i < r.size(); ++i)
          r[i] = 0;
        // set the diagonal elements
        r[0] = a[0];
        r[4] = a[1];
        r[8] = a[2];
        return vr;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::VECTOR_STRING) throw gError("FNDiagT::toC", "The diagMat() function may only be used with vector arguments");
      else
      {				
//         MSG_DEBUG("FNDiagT::toC", "writing");
        // the return type is always a tensor, independently of the argument 
        Variant vr(Variant::TENSOR_STRING);
        vector_string_t &r = vr.strings();
        vector_string_t a = vca.strings();
        
        // make sure that the off-diagonal elements are not undefined				
        for(size_t i = 0; i < r.size(); ++i)
          r[i] = string("(") + string("0") + string(")");
        // set the diagonal elements
        r[0] = "(" + a[0] + ")";
        r[4] = "(" + a[1] + ")";
        r[8] = "(" + a[2] + ")";
//         MSG_DEBUG("FNDiagT::toC", "vr[1] = " << vr[1]);
//         MSG_DEBUG("FNDiagT::toC", "vr[0] = " << vr[0]);
//         MSG_DEBUG("FNDiagT::toC", "vr[4] = " << vr[4]);
//         MSG_DEBUG("FNDiagT::toC", "vr[8] = " << vr[8]);
        return vr;							
      } 
    }									
};									
									
const FunctionNode_Register<FNDiagT> fn_diag_t(FunctionNode::UNARY_FUNCTION, 6, "diagMat");


class FNUnitV: public FNUnaryFunction	      
{							
  public:						
    FNUnitV(): FNUnaryFunction("idVec") {		
    }							

    virtual ~FNUnitV() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::SCALAR) throw gError("FNUnitV::value", "The idVec() function may only be used with scalar arguments");
      else
      {
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR);
        vector_double_t &vr = r.doubles();
        
        double d = vca.scalar();				
        for(size_t i = 0; i < vr.size(); ++i)
          vr[i] = d;
        return r;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::SCALAR_STRING) throw gError("FNUnitV::toC", "The idVec() function may only be used with scalar arguments");
      else
      {				
//         MSG_DEBUG("FNUnitV::toC", "writing");
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR_STRING);
        vector_string_t &vr = r.strings();
        
        for(size_t i = 0; i < vr.size(); ++i)
          vr[i] = "(" + vca.scalarString() + ")";
/*        for (size_t i = 0; i < vr.size(); ++i)
        MSG_DEBUG("FNUnitV::toC", "vr[" << i << "] = " << vr[i]);*/
        return r;							
      } 
    }									
};									
									
const FunctionNode_Register<FNUnitV> fn_unit_v(FunctionNode::UNARY_FUNCTION, 6, "idVec");


class FNIdMat: public FNUnaryFunction	      
{							
  public:						
    FNIdMat(): FNUnaryFunction("idMat") {		
    }							

    virtual ~FNIdMat() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::SCALAR) throw gError("FNIdMat::value", "The idMat() function may only be used with scalar arguments");
      else
      {
        // the return type is always a tensor, independently of the argument 
        Variant r(Variant::TENSOR);
        vector_double_t &vr = r.doubles();
        
        double d = vca.scalar();
        // make sure that the off-diagonal elements are not undefined				
        for(size_t i = 0; i < vr.size(); ++i)
          vr[i] = 0;
        // set the diagonal elements
        vr[0] = vr[4] = vr[8] = d;
        return r;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::SCALAR_STRING) throw gError("FNIdMat::toC", "The idMat() function may only be used with scalar arguments");
      else
      {				
//         MSG_DEBUG("FNIdMat::toC", "writing");
        // the return type is always a tensor, independently of the argument 
        Variant r(Variant::TENSOR_STRING);
        vector_string_t &vr = r.strings();
        
        // make sure that the off-diagonal elements are not undefined				
        for(size_t i = 0; i < vr.size(); ++i)
          vr[i] = string("(") + string("0") + string(")");
        // set the diagonal elements
        vr[0] = vr[4] = vr[8] = "(" + vca.scalarString() + ")";
//         MSG_DEBUG("FNIdMat::toC", "vr[1] = " << vr[1]);
//         MSG_DEBUG("FNIdMat::toC", "vr[0] = " << vr[0]);
//         MSG_DEBUG("FNIdMat::toC", "vr[4] = " << vr[4]);
//         MSG_DEBUG("FNIdMat::toC", "vr[8] = " << vr[8]);
        return r;							
      } 
    }									
};									
									
const FunctionNode_Register<FNIdMat> fn_id_mat(FunctionNode::UNARY_FUNCTION, 6, "idMat");


class FNQ: public FNUnaryFunction	      
{							
  public:						
    FNQ(): FNUnaryFunction("Q") {		
    }							

    virtual ~FNQ() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) {
        // the return type is always a scalar, independently from the argument 
        Variant r(Variant::SCALAR);
        // this one should have length one							
//         vector_double_t &vr = r.doubles();
        double& d = r.scalar();
//         assert(vr.size() == 1);
        /*vr[0]*/d = 0;
        vector_double_t va = vca.doubles();				
      
        for (size_t i = 0; i < va.size(); ++i) {				
          /*vr[0]*/d += va[i]*va[i];						
        }								
      
        return r;							
      } else								
        throw gError							
            ("FNQ::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	
          t == Variant::TENSOR_STRING) {				
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR_STRING);							
        // this one should have length one							
/*        vector_string_t &vr = r.strings();				
        assert(vr.size() == 1);*/
        /*vr[0]*/string d = "(";
        vector_string_t va = vca.strings();				
									
        for (size_t i = 0; i < va.size(); ++i) {				
          if(i) d/*vr[0]*/ += " + "; 
          d/*vr[0]*/ += "((" + va[i] + ")*(" + va[i] + "))";
        }								
        d/*vr[0]*/ += ")";
        r.scalarString() = d;
        					
        return r;							
          } else								
            throw gError							
                ("FNQ::toC",						
                 "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};									
									
const FunctionNode_Register<FNQ> fn_Q(FunctionNode::UNARY_FUNCTION, 6, "Q");



class FNStep: public FNUnaryFunction	      
{							
  public:						
    FNStep(): FNUnaryFunction("step") {		
    }							

    virtual ~FNStep() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) { 
        Variant r(t);							
        vector_double_t &vr = r.doubles();				
        vector_double_t va = vca.doubles();				
      
        for (size_t i = 0; i < va.size(); ++i) {				
          vr[i] = va[i] > 0 ? 1 : 0;						
        }								
      
        return r;							
      } else								
        throw gError							
            ("FNStep::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	
          t == Variant::TENSOR_STRING) {				
        Variant r(t);							
        vector_string_t &vr = r.strings();				
        vector_string_t va = vca.strings();				
									
        for (size_t i = 0; i < va.size(); ++i) {				
          vr[i] = "((" + va[i] + ") > 0 ? 1 : 0)";
        }								
									
        return r;							
          } else								
            throw gError							
                ("FNStep::toC",						
                 "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};									
									
const FunctionNode_Register<FNStep> fn_step(FunctionNode::UNARY_FUNCTION, 6, "step");

class FNStepVal: public FNUnaryFunction	      
{							
  public:						
    FNStepVal(): FNUnaryFunction("stpVal") {		
    }							

    virtual ~FNStepVal() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) { 
        Variant r(t);							
        vector_double_t &vr = r.doubles();				
        vector_double_t va = vca.doubles();				
      
        for (size_t i = 0; i < va.size(); ++i) {				
          vr[i] = va[i] > 0 ? va[i] : 0;						
        }								
      
        return r;							
      } else								
        throw gError							
            ("FNStepVal::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	
          t == Variant::TENSOR_STRING) {				
        Variant r(t);							
        vector_string_t &vr = r.strings();				
        vector_string_t va = vca.strings();				
									
        for (size_t i = 0; i < va.size(); ++i) {				
          vr[i] = "((" + va[i] + ") > 0 ? (" + va[i] + ") : 0)";
        }								
									
        return r;							
          } else								
            throw gError							
                ("FNStepVal::toC",						
                 "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};									
									
const FunctionNode_Register<FNStepVal> fn_step_val(FunctionNode::UNARY_FUNCTION, 6, "stpVal");


class FNTranspose: public FNUnaryFunction	      
{							
  public:						
    FNTranspose(): FNUnaryFunction("T") {		
    }							

    virtual ~FNTranspose() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      
      if (t != Variant::TENSOR) 
        throw gError("FNTranspose::value", "The transpose function T() may only be used with tensorial arguments.");
      else
      {
        // the return type is always a tensor 
        Variant vr(Variant::TENSOR);
        vector_double_t& r = vr.doubles();
        vector_double_t va = vca.doubles();				
        assert(va.size() == 9);
        r[0] = va[0]; r[1] = va[3]; r[2] = va[6];
        r[3] = va[1]; r[4] = va[4]; r[5] = va[7];
        r[6] = va[2]; r[7] = va[5]; r[8] = va[8];
        return vr;
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::TENSOR_STRING) {				
        Variant vr(Variant::TENSOR_STRING);							
        vector_string_t va = vca.strings();				
        vector_string_t& r = vr.strings();
        assert(va.size() == 9);
									
        r[0] = "(" + va[0] + ")"; r[1] = "(" + va[3] + ")"; r[2] = "(" + va[6] + ")";
        r[3] = "(" + va[1] + ")"; r[4] = "(" + va[4] + ")"; r[5] = "(" + va[7] + ")";
        r[6] = "(" + va[2] + ")"; r[7] = "(" + va[5] + ")"; r[8] = "(" + va[8] + ")";
                					
//         MSG_DEBUG("FNTranspose::toC", "r.size() = " << r.size());
        return vr;							
      } 
      else								
        throw gError							
            ("FNTranspose::toC",						
             "The transpose function T() may only be used with tensorial arguments");	
    }									
};									
									
const FunctionNode_Register<FNTranspose> fn_transpose(FunctionNode::UNARY_FUNCTION, 6, "T");


class FNTrace: public FNUnaryFunction	      
{							
  public:						
    FNTrace(): FNUnaryFunction("trace") {		
    }							

    virtual ~FNTrace() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      // the return type is always a scalar, independently of the argument 
      Variant r(Variant::SCALAR);
      
      if (t != Variant::TENSOR) throw gError("FNTrace::value", "The trace() function may only be used with tensorial arguments");
      else
      {
/*        // this one should have length one							
        vector_double_t &vr = r.doubles();
        assert(vr.size() == 1);*/
        double& d = r.scalar();
        vector_double_t va = vca.doubles();				
        assert(va.size() == 9);
        
        d/*vr[0]*/ = va[0]+va[4]+va[8];
        
      }
      return r;
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::TENSOR_STRING) {				
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR_STRING);							
        // this one should have length one							
/*        vector_string_t &vr = r.strings();				
        assert(vr.size() == 1);*/
        vector_string_t va = vca.strings();				
        assert(va.size() == 9);
									
        /*vr[0]*/r.scalarString() = "(" + va[0] + " + " + va[4] + " + " + va[8] + ")";
                					
                 return r;							
      } 
      else								
        throw gError							
            ("FNTrace::toC",						
             "The trace() function may only be used with tensorial arguments");	
    }									
};									
									
const FunctionNode_Register<FNTrace> fn_trace(FunctionNode::UNARY_FUNCTION, 6, "trace");


class FNUnitMat: public FNUnaryFunction	      
{							
  public:						
    FNUnitMat(): FNUnaryFunction("unitMat") {		
    }							

    virtual ~FNUnitMat() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::SCALAR) throw gError("FNUnitMat::value", "The unitMat() function may only be used with scalar arguments");
      else
      {
        // the return type is always a tensor, independently of the argument 
        Variant r(Variant::TENSOR);
        vector_double_t &vr = r.doubles();
        
        double d = vca.scalar();
        // make sure that the off-diagonal elements are not undefined				
        for(size_t i = 0; i < vr.size(); ++i)
          vr[i] = d;
        return r;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::SCALAR_STRING) throw gError("FNUnitMat::toC", "The unitMat() function may only be used with scalar arguments");
      else
      {				
//         MSG_DEBUG("FNUnitMat::toC", "writing");
        // the return type is always a tensor, independently of the argument 
        Variant r(Variant::TENSOR_STRING);
        vector_string_t &vr = r.strings();
        
        // make sure that the off-diagonal elements are not undefined				
        for(size_t i = 0; i < vr.size(); ++i)
          vr[i] = "(" + vca.scalarString() + ")";
        return r;							
      } 
    }									
};									
									
const FunctionNode_Register<FNUnitMat> fn_unit_mat(FunctionNode::UNARY_FUNCTION, 6, "unitMat");



// This is to provide a random number function for rns in [0,1]
class FNUran: public FNUnaryFunction	      
{							
  public:						
    FNUran(): FNUnaryFunction("uran") {		
    }							

    virtual ~FNUran() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR || t == Variant::VECTOR || t == Variant::TENSOR) { 
        Variant r(t);							
        vector_double_t &vr = r.doubles();				
        vector_double_t va = vca.doubles();				
      
        for (size_t i = 0; i < va.size(); ++i) {				
          vr[i] = (double)rand()/(double)RAND_MAX;
        }								
      
        return r;							
      } else								
        throw gError							
            ("FNStep::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING || t == Variant::VECTOR_STRING ||	
          t == Variant::TENSOR_STRING) {				
        Variant r(t);							
        vector_string_t &vr = r.strings();				
        vector_string_t va = vca.strings();				
									
        for (size_t i = 0; i < va.size(); ++i) {				
          vr[i] = "((double)rand()/(double)RAND_MAX)";
        }								
									
        return r;							
          } else								
            throw gError							
                ("FNStep::toC",						
                 "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};	

const FunctionNode_Register<FNUran> fn_uran(FunctionNode::UNARY_FUNCTION, 6, "uran");


class FNUnitVX: public FNUnaryFunction	      
{							
  public:						
    FNUnitVX(): FNUnaryFunction("uVecX") {		
    }							

    virtual ~FNUnitVX() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::SCALAR) throw gError("FNUnitVX::value", "The uVecX() function may only be used with scalar arguments");
      else
      {
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR);
        vector_double_t &vr = r.doubles();
        
        double d = vca.scalar();				
        vr[0] = d;
        vr[1] = 0;
        vr[2] = 0;
        return r;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::SCALAR_STRING) throw gError("FNUnitVX::toC", "The uVecX() function may only be used with scalar arguments");
      else
      {				
//         MSG_DEBUG("FNUnitV::toC", "writing");
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR_STRING);
        vector_string_t &vr = r.strings();
        
        vr[0] = "(" + vca.scalarString() + ")";
        vr[1] = "(0)";
        vr[2] = "(0)";
        return r;							
      } 
    }									
};									
									
const FunctionNode_Register<FNUnitVX> fn_uvec_x(FunctionNode::UNARY_FUNCTION, 6, "uVecX");

class FNUnitVY: public FNUnaryFunction	      
{							
  public:						
    FNUnitVY(): FNUnaryFunction("uVecY") {		
    }							

    virtual ~FNUnitVY() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::SCALAR) throw gError("FNUnitVY::value", "The uVecY() function may only be used with scalar arguments");
      else
      {
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR);
        vector_double_t &vr = r.doubles();
        
        double d = vca.scalar();				
        vr[1] = d;
        vr[2] = 0;
        vr[0] = 0;
        return r;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::SCALAR_STRING) throw gError("FNUnitVY::toC", "The uVecY() function may only be used with scalar arguments.");
      else
      {				
//         MSG_DEBUG("FNUnitV::toC", "writing");
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR_STRING);
        vector_string_t &vr = r.strings();
        
        vr[1] = "(" + vca.scalarString() + ")";
        vr[2] = "(0)";
        vr[0] = "(0)";
        return r;							
      } 
    }									
};									
									
const FunctionNode_Register<FNUnitVY> fn_uvec_y(FunctionNode::UNARY_FUNCTION, 6, "uVecY");

class FNUnitVZ: public FNUnaryFunction	      
{							
  public:						
    FNUnitVZ(): FNUnaryFunction("uVecZ") {		
    }							

    virtual ~FNUnitVZ() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t != Variant::SCALAR) throw gError("FNUnitVZ::value", "The uVecZ() function may only be used with scalar arguments");
      else
      {
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR);
        vector_double_t &vr = r.doubles();
        
        double d = vca.scalar();				
        vr[2] = d;
        vr[0] = 0;
        vr[1] = 0;
        return r;							
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t != Variant::SCALAR_STRING) throw gError("FNUnitVZ::toC", "The uVecZ() function may only be used with scalar arguments");
      else
      {				
//         MSG_DEBUG("FNUnitV::toC", "writing");
        // the return type is always a vector, independently of the argument 
        Variant r(Variant::VECTOR_STRING);
        vector_string_t &vr = r.strings();
        
        vr[2] = "(" + vca.scalarString() + ")";
        vr[0] = "(0)";
        vr[1] = "(0)";
        return r;							
      } 
    }									
};									
									
const FunctionNode_Register<FNUnitVZ> fn_uvec_z(FunctionNode::UNARY_FUNCTION, 6, "uVecZ");


class FNxyMat: public FNUnaryFunction	      
{							
  public:						
    FNxyMat(): FNUnaryFunction("xyMat") {		
    }							

    virtual ~FNxyMat() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      
      if (t != Variant::TENSOR) 
        throw gError("FNxyMat::value", "The function xyMat() may only be used with tensorial arguments.");
      else
      {
        // the return type is always a tensor 
        Variant vr(Variant::TENSOR);
        vector_double_t& r = vr.doubles();
        vector_double_t va = vca.doubles();				
        assert(va.size() == 9);
        r[0] = va[0]; r[1] = va[1]; r[2] = 0;
        r[3] = va[3]; r[4] = va[4]; r[5] = 0;
        r[6] = 0; r[7] = 0; r[8] = 0;
        return vr;
      }
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::TENSOR_STRING) {				
        Variant vr(Variant::TENSOR_STRING);							
        vector_string_t va = vca.strings();				
        vector_string_t& r = vr.strings();
        assert(va.size() == 9);
									
        r[0] = "(" + va[0] + ")"; r[1] = "(" + va[1] + ")"; r[2] = "(" + string("0") + ")";
        r[3] = "(" + va[3] + ")"; r[4] = "(" + va[4] + ")"; r[5] = "(" + string("0") + ")";
        r[6] = "(" + string("0") + ")"; r[7] = "(" + string("0") + ")"; r[8] = "(" + string("0") + ")";
                					
//         MSG_DEBUG("FNxyMat::toC", "r.size() = " << r.size());
        return vr;							
      } 
      else								
        throw gError							
            ("FNxyMat::toC",						
             "The function xyMat() may only be used with tensorial arguments");	
    }									
};									
									
const FunctionNode_Register<FNxyMat> fn_xy_mat(FunctionNode::UNARY_FUNCTION, 6, "xyMat");



class FNXPos: public FNUnaryFunction	      
{							
  public:						
    FNXPos(): FNUnaryFunction("xCoord") {		
    }							

    virtual ~FNXPos() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR) throw gError("FNXPos::value", "The xCoord() function may not be used with scalar arguments");
      
      if (t == Variant::TENSOR) throw gError("FNXPos::value", "The xCoord() function may not be used with tensor arguments");
      
      if (Variant::VECTOR) {
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR);
        // this one should have length one							
        vector_double_t &vr = r.doubles();
        assert(vr.size() == 1);
        vector_double_t va = vca.doubles();				
        assert(va.size() == 3);
        
        vr[0] = va[0];
        
        return r;							
      } else								
        throw gError							
            ("FNXPos::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING) throw gError("FNXPos::value", "The xCoord() function may not be used with scalar arguments");
      
      if (t == Variant::TENSOR_STRING) throw gError("FNXPos::value", "The xCoord() function may not be used with tensor arguments");
      
      if (t == Variant::VECTOR_STRING) {				
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR_STRING);							
        // this one should have length one							
        vector_string_t &vr = r.strings();				
        assert(vr.size() == 1);
        vector_string_t va = vca.strings();				
        assert(va.size() == 3);
									
        vr[0] = "(" + va[0] + ")";
                					
        return r;							
      } 
      else								
        throw gError							
            ("FNXPos::toC",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};									
									
const FunctionNode_Register<FNXPos> fn_xpos(FunctionNode::UNARY_FUNCTION, 6, "xCoord");



class FNYPos: public FNUnaryFunction	      
{							
  public:						
    FNYPos(): FNUnaryFunction("yCoord") {		
    }							

    virtual ~FNYPos() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR) throw gError("FNYPos::value", "The yCoord() function may not be used with scalar arguments");
      
      if (t == Variant::TENSOR) throw gError("FNYPos::value", "The yCoord() function may not be used with tensor arguments");
      
      if (Variant::VECTOR) {
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR);
        // this one should have length one							
        vector_double_t &vr = r.doubles();
        assert(vr.size() == 1);
        vector_double_t va = vca.doubles();				
        assert(va.size() == 3);
        
        vr[0] = va[1];
        
        return r;							
      } else								
        throw gError							
            ("FNYPos::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING) throw gError("FNYPos::value", "The yCoord() function may not be used with scalar arguments");
      
      if (t == Variant::TENSOR_STRING) throw gError("FNYPos::value", "The yCoord() function may not be used with tensor arguments");
      
      if (t == Variant::VECTOR_STRING) {				
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR_STRING);							
        // this one should have length one							
        vector_string_t &vr = r.strings();				
        assert(vr.size() == 1);
        vector_string_t va = vca.strings();				
        assert(va.size() == 3);
									
        vr[0] = "(" + va[1] + ")";
                					
        return r;							
      } 
      else								
        throw gError							
            ("FNYPos::toC",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};									
									
const FunctionNode_Register<FNYPos> fn_ypos(FunctionNode::UNARY_FUNCTION, 6, "yCoord");



class FNZPos: public FNUnaryFunction	      
{							
  public:						
    FNZPos(): FNUnaryFunction("zCoord") {		
    }							

    virtual ~FNZPos() {				
    }							
      							
    virtual Variant value() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->value();					
									
      t = vca.typeId();							
									
      if (t == Variant::SCALAR) throw gError("FNZPos::value", "The zCoord() function may not be used with scalar arguments");
      
      if (t == Variant::TENSOR) throw gError("FNZPos::value", "The zCoord() function may not be used with tensor arguments");
      
      if (Variant::VECTOR) {
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR);
        // this one should have length one							
        vector_double_t &vr = r.doubles();
        assert(vr.size() == 1);
        vector_double_t va = vca.doubles();				
        assert(va.size() == 3);
        
        vr[0] = va[2];
        
        return r;							
      } else								
        throw gError							
            ("FNZPos::value",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
    									
    virtual Variant toC() const {					
      Variant::variant_type_t t;					
      Variant vca = m_a->toC();						
									
      t = vca.typeId();							
    
      if (t == Variant::SCALAR_STRING) throw gError("FNZPos::value", "The zCoord() function may not be used with scalar arguments");
      
      if (t == Variant::TENSOR_STRING) throw gError("FNZPos::value", "The zCoord() function may not be used with tensor arguments");
      
      if (t == Variant::VECTOR_STRING) {				
        // the return type is always a scalar, independently of the argument 
        Variant r(Variant::SCALAR_STRING);							
        // this one should have length one							
        vector_string_t &vr = r.strings();				
        assert(vr.size() == 1);
        vector_string_t va = vca.strings();				
        assert(va.size() == 3);
									
        vr[0] = "(" + va[2] + ")";
                					
        return r;							
      } 
      else								
        throw gError							
            ("FNZPos::toC",						
             "Can't handle Variant type id " + ObjToString(t) + ".");	
    }									
};									
									
const FunctionNode_Register<FNZPos> fn_zpos(FunctionNode::UNARY_FUNCTION, 6, "zCoord");
