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


#ifndef __PCA_MATRIX_INVERSE_H
#define __PCA_MATRIX_INVERSE_H 

#include <gsl/gsl_math.h>

#include "gsl_helper.h"

#include "particle_cache.h"


/*!
 * Computes the inverse of a given symmetric square matrix 
 * The members inherited from class Symbol (m_symbolName, m_offset, m_datatype) 
 * will be used for the computed inverse.
 */
class PCaMatrixInverse: public ParticleCache
{
  protected:
  
  /*!
   * Tag offset of the given tensor
   */
    size_t m_tensor_offset;

    /*!
     * Symbol of the given tensor
     */
    string m_tensor_symbol;
    
    /*!
     * Used for copying the given matrix. This is OK/necessary because 
     * gsl_linalg_LU_decomp modifies it during LU decomposition
     */
    gsl_matrix* m_working_mat;
    
    /*!
     * Used for storing the inverse. We do not allocate memory for the pointer double*
     * because it will be directed to the right adress in the particle tag. The other 
     * information is set by hand in the constructor.
     */
    gsl_matrix m_inverse;
    
    /*!
     * The workspace for the computation of the LU decomposition
     */
    gsl_permutation *m_permutation;
    
    /*!
     * Initialise the property list
     */
    virtual void init();

    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new PCaMatrixInverse(*this);
    }
    
    /*!
     * Adds the expressions used by this \a Symbol to the given list. 
     * @param usedSymbols List to be filled with own instances of \a TypedValue
     */
    virtual void addMyUsedSymbolsTo(typed_value_list_t& usedSymbols)
    {

    }

    /*!
     * Returns the strings of those \a Symbols that the given class depends on
     * due to hard-coded reasons (not due to runtime compiled expressions).
     * @param usedSymbols List to add the strings to.
     */
    virtual void addMyHardCodedDependenciesTo(list<string>& usedSymbols) const
    {
      usedSymbols.push_back(m_tensor_symbol);
    }



 public:
    
    /*!
   * Constructor for node hierarchy
     */
    PCaMatrixInverse(/*Node*/Simulation* parent);
 
  /*!
     * Destructor
   */
    virtual ~PCaMatrixInverse();

  /*!
     * Compute the inverse
   */
    virtual void computeCacheFor(Particle* p) {/*
      p->tag.doubleByOffset(m_offset) += m_wf->interpolate(NULL, p->r) / p->tag.doubleByOffset(m_density_offset);*/
      
//         MSG_DEBUG("PCaMatrixInverse::computeCacheFor", "shear tensor: " << p->tag.tensorByOffset(8));
      
//         MSG_DEBUG("PCaMatrixInverse::computeCacheFor", "original tensor: " << p->tag.tensorByOffset(m_tensor_offset));
      
      
      // assign the right memory adress to m_inverse
      m_inverse.data = ((p->tag.tensorByOffset(m_offset)).tensor);

      
      // load the tensor
      tensor_t& tensor = p->tag.tensorByOffset(m_tensor_offset);
//       size_t k = 0;
//       MSG_DEBUG("PCaMatrixInverse::computeCacheFor", "copied reference of tensor: " << tensor);
      
//       MSG_DEBUG("PCaMatrixInverse::computeCacheFor", "GSL-working-tensor: ");
      for(size_t i = 0; i < SPACE_DIMS; ++i)
      {
        for(size_t j = 0; j < SPACE_DIMS; ++j)
        {
          GSL_MATRIX_SET(m_working_mat, i, j, tensor(i,j));
        
//           cout << GSL_MATRIX_GET(m_working_mat, i, j) << " ";
        }
//         cout << endl;
      } 
            
      // Compute the LU decomposition
      int signum;
      int error = gsl_linalg_LU_decomp(m_working_mat, m_permutation, &signum);
      if(error)
/*         throw gError("PCaMatrixInverse::computeCacheFor", "gsl_linalg_LU_decomp failed with error " + string(gsl_strerror(error))); */
        throw gError("PCaMatrixInverse::computeCacheFor", "gsl_linalg_LU_decomp failed with return value " + ObjToString(error));
      
      // Compute the inverse
      error = gsl_linalg_LU_invert(m_working_mat, m_permutation, &m_inverse);
      if(error) {
	MSG_DEBUG("PCaMatrixInverse::computeCacheFor", "ERROR for Marix: " << tensor << ", particle: " << p->mySlot << ", colour: " << p->c << ":");
	if(gsl_linalg_LU_det(m_working_mat, signum) == 0.0)
	  throw gError("PCaMatrixInverse::computeCacheFor", "failed because matrix is singular.");
	else
	  throw gError("PCaMatrixInverse::computeCacheFor", "gsl_linalg_LU_invert failed with error " + string(gsl_strerror(error))); 
/*         throw gError("PCaMatrixInverse::computeCacheFor", "gsl_linalg_LU_invert failed with return value " + ObjToString(error)); */

      }

    }

  /*!
     * Take steps necessary to register this calculator
   */
    virtual void registerWithParticle() {}

  /*!
     * Does this calculator equal \a c?
     * @param c Other calculator
   */
    virtual bool operator==(const ParticleCache &c) const {
      if (typeid(c) == typeid(*this)) {
        PCaMatrixInverse *cc = (PCaMatrixInverse*) &c;

        return
            m_offset == cc->m_offset && m_tensor_offset == cc->m_tensor_offset &&  m_tensor_symbol == cc->m_tensor_symbol && m_stage == cc->m_stage && m_offset == cc->m_offset && m_colour == cc->m_colour;
      } else
        return false;
    }
      
    /*!
     * Return the name of the computed symbols to be used in other expressions.
     */
    virtual list<string> mySymbolNames()
    {
      list<string> temp;
      assert(temp.empty());
      temp.push_back(m_symbolName);
      return temp;
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();
        
};

#endif
