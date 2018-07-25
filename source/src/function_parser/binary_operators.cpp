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


#include "binary_operators.h"

#include "function_parser.h"
#include "typed_value.h"

/* BinaryOperator */

FNBinaryOperator::FNBinaryOperator(string name): FunctionNode(name), m_a(NULL), m_b(NULL)
{
}


FNBinaryOperator::~FNBinaryOperator()
{
  if (!dynamic_cast<TypedValue*>(m_a))
    delete m_a;
  if (!dynamic_cast<TypedValue*>(m_b))
    delete m_b;
}


void FNBinaryOperator::setBinary(FunctionNode *a, FunctionNode *b)
{
  if (!dynamic_cast<TypedValue*>(m_a))
    delete m_a;
  if (!dynamic_cast<TypedValue*>(m_b))
    delete m_b;

  m_a = a;
  m_b = b;
}



/*!
 * The addition operator. Acts on scalars, vectors, and tensors.
 */
class FNAddition: public FNBinaryOperator {
public:
  /*!
   * Constructor
   */
  FNAddition(): FNBinaryOperator("+") {
  }

  /*!
   * Destructor
   */
  virtual ~FNAddition() {
  }

  /*!
   * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
  virtual Variant value() const {
    FOR_EACH_DOUBLE2
      (FNAddition::value,
       m_a->value(),
       m_b->value(),
       vr[i] = va[i] + vb[i];
       );
  }

  /*!
   * Calculates the expression in C notation.
   */
  virtual Variant toC() const {
    FOR_EACH_STRING2
      (FNAddition::toC,
       m_a->toC(),
       m_b->toC(),
       vr[i] = "(" + va[i] + "+" + vb[i] + ")";
       );
  }
};

const FunctionNode_Register<FNAddition> fn_addition(FunctionNode::BINARY_OPERATOR, 1, "+");



/*!
 * The subtraction operator. Acts on scalar, vectors, and tensors.
 */
class FNSubtraction: public FNBinaryOperator {
public:
  /*!
   * Constructor
   */
  FNSubtraction(): FNBinaryOperator("-") {
  }

  /*!
   * Destructor
   */
  virtual ~FNSubtraction() {
  }

  /*!
   * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
  virtual Variant value() const {
    FOR_EACH_DOUBLE2
      (FNSubstraction::value,
       m_a->value(),
       m_b->value(),
       
       vr[i] = va[i] - vb[i];
       );
  }

  /*!
   * Calculates the expression in C notation.
   */
  virtual Variant toC() const {
    FOR_EACH_STRING2
      (FNAddition::toC,
       m_a->toC(),
       m_b->toC(),
       vr[i] = "(" + va[i] + "-" + vb[i] + ")";
       );
  }
};

const FunctionNode_Register<FNSubtraction> fn_subtraction(FunctionNode::BINARY_OPERATOR, 1, "-");



/*!
 * The multiplication operator. Acts on scalars, vectors, and tensors.
 * For vectors and tensors the elements are multiplied component wise.
 * An inner product is given by the contraction ':' operator.
 */
class FNMultiplication: public FNBinaryOperator { 
public:
  /*!
   * Constructor
   */
  FNMultiplication(): FNBinaryOperator("*"){
  }

  /*!
   * Destructor
   */
  virtual ~FNMultiplication() {
  }

  /*!
   * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
  virtual Variant value() const {
    Variant vara = m_a->value();
    Variant varb = m_b->value();

    if (vara.typeId() == Variant::SCALAR) {
      FOR_EACH_DOUBLE
	(FNMultiplication::value,
	 varb,
	 vr[i] = vara.scalar()*va[i]
	 );
    } else if (varb.typeId() == Variant::SCALAR) {
      FOR_EACH_DOUBLE
	(FNMultiplication::value,
	 vara,
	 vr[i] = varb.scalar()*va[i]
	 );
    } else {
      FOR_EACH_DOUBLE2
	(FNMultiplication::value,
	 vara,
	 varb,
	 
	 vr[i] = va[i] * vb[i];
	 );
    }
  }

  /*!
   * Calculates the expression in C notation.
   */
  virtual Variant toC() const {
    Variant vara = m_a->toC();
    Variant varb = m_b->toC();

    if (vara.typeId() == Variant::SCALAR_STRING) {
      FOR_EACH_STRING
	(FNMultiplication::toC,
	 varb,
	 vr[i] = "(" + vara.scalarString() + "*" + va[i] + ")";
	 );
    } else if (varb.typeId() == Variant::SCALAR_STRING) {
      FOR_EACH_STRING
	(FNMultiplication::toC,
	 vara,
	 vr[i] = "(" + varb.scalarString() + "*" + va[i] + ")";
	 );
    } else {
      FOR_EACH_STRING2
	(FNMultiplication::toC,
	 m_a->toC(),
	 m_b->toC(),
	 vr[i] = "(" + va[i] + "*" + vb[i] + ")";
	 );
    }
  }
};

const FunctionNode_Register<FNMultiplication> fn_multiplication(FunctionNode::BINARY_OPERATOR, 2, "*");



/*!
 * The division operator. Acts on scalars, vectors, and tensors.
 * For vectors and tensors the elements are divided component wise.
 */
class FNDivision: public FNBinaryOperator {
public:
  /*!
   * Constructor
   */
  FNDivision(): FNBinaryOperator("/") {
  }
  
  /*!
   * Destructor
   */
  virtual ~FNDivision() {
  }

  /*!
   * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
  virtual Variant value() const {
    Variant vara = m_a->value();
    Variant varb = m_b->value();

    if (varb.typeId() == Variant::SCALAR) {
      FOR_EACH_DOUBLE
	(FNDivision::value,
	 vara,
	 vr[i] = va[i]/varb.scalar()
	 );
    } else {
      FOR_EACH_DOUBLE2
	(FNDivision::value,
	 vara,
	 varb,
	 
	 vr[i] = va[i] / vb[i];
	 );
    }
  }

  /*!
   * Calculates the expression in C notation.
   */
  virtual Variant toC() const {
    Variant vara = m_a->toC();
    Variant varb = m_b->toC();

    if (varb.typeId() == Variant::SCALAR_STRING) {
      FOR_EACH_STRING
	(FNDivision::toC,
	 vara,
	 vr[i] = "(" + va[i] + "/" + varb.scalarString() + ")";	 
	 );
    } else {
      FOR_EACH_STRING2
	(FNDivision::toC,
	 m_a->toC(),
	 m_b->toC(),
	 vr[i] = "(" + va[i] + "/" + vb[i] + ")";
	 );
    }
  }
};

const FunctionNode_Register<FNDivision> fn_division(FunctionNode::BINARY_OPERATOR, 3, "/");



/*!
 * The contraction operator. For vectors and tensors all elements are contracted,
 * i.e. multiplied and summed. This gives, e.g., the dot product for two vectors.
 */
class FNContraction: public FNBinaryOperator {
  public:
  /*!
   * Constructor
   */
    FNContraction(): FNBinaryOperator(":") {
    }

  /*!
     * Destructor
   */
    virtual ~FNContraction() {
    }

  /*!
     * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
    virtual Variant value() const {
      Variant va = m_a->value();
      Variant vb = m_b->value();

      if ((va.typeId() == Variant::SCALAR && vb.typeId() == Variant::SCALAR) || (va.typeId() == Variant::VECTOR && vb.typeId() == Variant::VECTOR) ||
           (va.typeId() == Variant::TENSOR && vb.typeId() == Variant::TENSOR)) {
        Variant r(Variant::SCALAR);
        double d = 0;
        vector_double_t a = va.doubles();
        vector_double_t b = vb.doubles();
      
        for (size_t i = 0; i < a.size(); ++i)
          d += a[i]*b[i];
      
        r.scalar() = d;

        return r;
           }
           else if(va.typeId() == Variant::TENSOR && vb.typeId() == Variant::VECTOR)
           {
             Variant vr(Variant::VECTOR);
             vector_double_t a = va.doubles();
             vector_double_t b = vb.doubles();
             vector_double_t& r = vr.doubles();
             size_t dim = b.size();
             for (size_t i = 0; i < dim; ++i)
               r[i] = 0;
             for (size_t i = 0; i < dim; ++i)
             {
               for (size_t j = 0; j < dim; ++j)
               {
                 r[j] += a[i+j*dim]*b[i];
               }      
             }
             return vr;
           }
           else 
             throw gError
                 ("FNContraction::value",
                  "Incorrect types for operand ':'. Either, both operands must be of the same type or the first operand must be a tensor and the second one a vector.");
    }

  /*!
     * Calculates the expression in C notation.
   */
    virtual Variant toC() const {
      Variant va = m_a->toC();
      Variant vb = m_b->toC();
    
      if ((va.typeId() == Variant::VECTOR_STRING && vb.typeId() == Variant::VECTOR_STRING) ||
           (va.typeId() == Variant::TENSOR_STRING && vb.typeId() == Variant::TENSOR_STRING)) {
        Variant r(Variant::SCALAR_STRING);
        string d = "";
        vector_string_t a = va.strings();
        vector_string_t b = vb.strings();
      
        for (size_t i = 0; i < a.size(); ++i) 
        {
          if (i) d += "+";
          d += a[i] + "*" + b[i];
        }
      
        r.scalarString() = "(" + d + ")";
      
        return r;
           }
           else if(va.typeId() == Variant::TENSOR_STRING && vb.typeId() == Variant::VECTOR_STRING)
           {
	     
             Variant r(Variant::VECTOR_STRING);
             vector_string_t& vr = r.strings();
	     vector_string_t a = va.strings();
	     vector_string_t b = vb.strings();
	     
	     size_t dim = b.size();
	     for (size_t i = 0; i < dim; ++i)
	       vr[i] = "(";
	     for (size_t i = 0; i < dim; ++i)
	       {
		 for (size_t j = 0; j < dim; ++j)
		   {
		     if(i) vr[j] += "+";
		     vr[j] += a[i+j*dim] + "*" + b[i];
		   }      
	       }
	     for (size_t i = 0; i < dim; ++i)
	       vr[i] += ")";
	     return r;
           }
           else 
             throw gError
                 ("FNContraction::toC",
                  "Incorrect types for operand ':'. Either, both operands must be of the same type or the first operand must be a tensor and the second one a vector.");
    }
};

const FunctionNode_Register<FNContraction> fn_contraction(FunctionNode::BINARY_OPERATOR, 3, ":");

/*!
 * The "dot" operator.
 */
class FNDot: public FNBinaryOperator {
  public:
  /*!
   * Constructor
   */
    FNDot(): FNBinaryOperator("°") {
    }

  /*!
     * Destructor
   */
    virtual ~FNDot() {
    }

  /*!
     * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
    virtual Variant value() const {
      Variant va = m_a->value();
      Variant vb = m_b->value();

      if ((va.typeId() == Variant::TENSOR && vb.typeId() == Variant::TENSOR)) 
      {
        Variant vr(Variant::TENSOR);
        vector_double_t& r = vr.doubles();
        vector_double_t a = va.doubles();
        vector_double_t b = vb.doubles();
        for (size_t i = 0; i < SPACE_DIMS*SPACE_DIMS; ++i)
          r[i] = 0;
      
        for (size_t i = 0; i < SPACE_DIMS; ++i)
          for (size_t j = 0; j < SPACE_DIMS; ++j)
            for (size_t k = 0; k < SPACE_DIMS; ++k)
              r[j+i*SPACE_DIMS] += a[k+i*SPACE_DIMS]*b[j+k*SPACE_DIMS];

        return vr;
      }
      else 
        throw gError
            ("FNDot::value",
             "Incorrect types for operand '°'. Currently, only two matrices are allowed.");
    }

  /*!
     * Calculates the expression in C notation.
   */
    virtual Variant toC() const {
      Variant va = m_a->toC();
      Variant vb = m_b->toC();
    
      if ((va.typeId() == Variant::TENSOR_STRING && vb.typeId() == Variant::TENSOR_STRING)) {
        Variant vr(Variant::TENSOR_STRING);
        vector_string_t& r = vr.strings();
        vector_string_t a = va.strings();
        vector_string_t b = vb.strings();
      
        for (size_t i = 0; i < SPACE_DIMS*SPACE_DIMS; ++i)
          r[i] = "("; 
        
        for (size_t i = 0; i < SPACE_DIMS; ++i) 
          for (size_t j = 0; j < SPACE_DIMS; ++j) 
            for (size_t k = 0; k < SPACE_DIMS; ++k) 
        {
          if(k == 0) r[j+i*SPACE_DIMS] += "(";
          if(k) r[j+i*SPACE_DIMS] += " + ";
          r[j+i*SPACE_DIMS] += "(" + a[k+i*SPACE_DIMS] + "*" + b[j+k*SPACE_DIMS] + ")";
          if(k == SPACE_DIMS-1) r[j+i*SPACE_DIMS] += ")";
        }
      
        for (size_t i = 0; i < SPACE_DIMS*SPACE_DIMS; ++i)
          r[i] += ")"; 
      
        return vr;
      }
      else 
        throw gError
            ("FNDot::toC",
             "Incorrect types for operand '°'. Currently, only two matrices are allowed.");
    }
};

const FunctionNode_Register<FNDot> fn_dot(FunctionNode::BINARY_OPERATOR, 3, "°");


/*!
 * The "outer-product" operator.
 */
class FNOuter: public FNBinaryOperator {
  public:
  /*!
   * Constructor
   */
    FNOuter(): FNBinaryOperator("@") {
    }

  /*!
     * Destructor
   */
    virtual ~FNOuter() {
    }

  /*!
     * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
    virtual Variant value() const {
      Variant va = m_a->value();
      Variant vb = m_b->value();

      if ((va.typeId() == Variant::VECTOR && vb.typeId() == Variant::VECTOR)) 
      {
        Variant vr(Variant::TENSOR);
        vector_double_t& r = vr.doubles();
        vector_double_t a = va.doubles();
        vector_double_t b = vb.doubles();
        for (size_t i = 0; i < SPACE_DIMS; ++i)
          for (size_t j = 0; j < SPACE_DIMS; ++j)
              r[j+i*SPACE_DIMS] = a[i]*b[j];

        return vr;
      }
      else 
        throw gError
            ("FNOuter::value",
             "Incorrect types for operand '@'. Only two vectors are allowed.");
    }

  /*!
     * Calculates the expression in C notation.
   */
    virtual Variant toC() const {
      Variant va = m_a->toC();
      Variant vb = m_b->toC();
    
      if ((va.typeId() == Variant::VECTOR_STRING && vb.typeId() == Variant::VECTOR_STRING)) {
        Variant vr(Variant::TENSOR_STRING);
        vector_string_t& r = vr.strings();
        vector_string_t a = va.strings();
        vector_string_t b = vb.strings();
      
        for (size_t i = 0; i < SPACE_DIMS*SPACE_DIMS; ++i)
          r[i] = "("; 
        
        for (size_t i = 0; i < SPACE_DIMS; ++i) 
          for (size_t j = 0; j < SPACE_DIMS; ++j) 
	    {
	      r[j+i*SPACE_DIMS] += a[i] + "*" + b[j];
	    }
      
        for (size_t i = 0; i < SPACE_DIMS*SPACE_DIMS; ++i)
          r[i] += ")"; 
      
        return vr;
      }
      else 
        throw gError
            ("FNouter::toC",
             "Incorrect types for operand '@'. Only two vectors are allowed.");
    }
};

const FunctionNode_Register<FNOuter> fn_outer(FunctionNode::BINARY_OPERATOR, 3, "@");

/*!
 * The power operator.
 * Currently it is only defined for two scalars. Is there a definition for higher
 * dimensional quantities?
 */
class FNPower: public FNBinaryOperator {
public:
  /*!
   * Constructor
   */
  FNPower(): FNBinaryOperator("^") {
  }

  /*!
   * Destructor
   */
  virtual ~FNPower() {
  }

  /*!
   * Calculates the *value* of this operation, i.e. evaluate the expression.
   */
  virtual Variant value() const {
    Variant va = m_a->value();
    Variant vb = m_b->value();

    if (va.typeId() == Variant::SCALAR && vb.typeId() == Variant::SCALAR) {
      Variant r(Variant::SCALAR);

      r.scalar() = pow(va.scalar(), vb.scalar());

      return r;
    } else 
      throw gError
	("FNPower::value",
	 "Both operands must be scalars for the power operation.");
  }

  /*!
   * Calculates the expression in C notation.
   */
  virtual Variant toC() const {

    Variant va = m_a->toC();
    Variant vb = m_b->toC();

    if (va.typeId() == Variant::SCALAR_STRING && vb.typeId() == Variant::SCALAR_STRING) {
      
      Variant r(Variant::SCALAR_STRING);

      // This message could be obsolete due to the added (2018-02-12)
      // try .. catch that follow
      // FUNCTION_PARSER_LOG("FNPower::toC()", "Now I check the second argument of the '^'-operator. If I boil out with s.th. like a seg.-fault, it is neither a real number nor an integer and you should check your expression.");    
      
      // if m_b is a non-primitive scalar, i.e., an expression with
      // Symbols, then m_b->value() will throw an exception. This one
      // we try to catch and see if we can create a "pow(..)" expression
      // nonetheless
      try {
      
      Variant vValB = m_b->value();

      int tmpInt = int(vValB.scalar());

      // is the second argument an integer?
      if(tmpInt == vValB.scalar())
      {
	// yes, so we can do it the fast way
        if(tmpInt > 0)
        {
          r.scalarString() = "(";
          r.scalarString() += va.scalarString();
          for(int i = 1; i < tmpInt; ++i) 
            r.scalarString() += "*" + va.scalarString();
          r.scalarString() += ")";
        }
        else
        {
	  // exponent is zero, so result is "1"
          if(tmpInt == 0) r.scalarString() = "(1)";
          else
          {
            // exponent b is negative, so compute 1/(a^(-b))
            tmpInt = -tmpInt;
            r.scalarString() = "(1.0/(" + va.scalarString();
            for(int i = 1; i < tmpInt; ++i)
              r.scalarString() += "*" + va.scalarString();
            r.scalarString() += "))";
          }
        }
      }
      else
        // unfortunately the exponent is not an integer, so we have to take the slow 
        // pow() function
        r.scalarString() = "(pow(" + va.scalarString() + ", " + vb.scalarString() + "))";
      
      return r;

      } // end of try{..}
      catch(gError &err) {

	FUNCTION_PARSER_LOG("FNPower::toC()", "Now I try to make a function with variables out of the second argument of the '^'-operator. If I boil out with s.th. like a seg.-fault, then I am probably not yet clever enough to digest your expression '" << vb.scalarString() << "'.");    

	// try if we get something useful with the pow function
        r.scalarString() = "(pow(" + va.scalarString() + ", " + vb.scalarString() + "))";

	return r;
      }

      } else 
      throw gError
	("FNPower::toC",
	 "Both operands must be scalars for the power operand '^'.");
  }
};

const FunctionNode_Register<FNPower> fn_power(FunctionNode::BINARY_OPERATOR, 4, "^");
