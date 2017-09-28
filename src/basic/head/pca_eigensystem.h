/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2017, 
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


#ifndef __PCA_EIGENSYSTEM_H
#define __PCA_EIGENSYSTEM_H 

#include "particle_cache.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// the following definitions are done because I do not know how to turn range-checking 
// off for the GSL
#define GSL_VECTOR_SET(v, i, x) (v)->data[i*(v)->stride] = x
#define GSL_VECTOR_GET(v, i) (v)->data[i*(v)->stride]
#define GSL_MATRIX_SET(m, i, j, x) (m)->data[i * (m)->tda + j] = x
#define GSL_MATRIX_GET(m, i, j) (m)->data[i * (m)->tda + j]

/*!
 * Computes the eigensystem, i.e., the eigenvalues and corresponding 
 * eigenvectors for a given symmetric square matrix 
 * The members inherited from class Symbol (m_symbolName, m_offset, m_datatype) 
 * will be used for the vector of eigenvalues. For the matrix of eigenvectors, 
 * additional members are added
 */
class PCaEigensystem: public ParticleCache
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
     * Tag offset of matrix of eigenvectors
   */
    size_t m_evecs_offset;

    /*!
     * Symbol of the given tensor
     */
    string m_evecs_symbol;
    
    /*!
     * Used for copying the current matrix. This is OK/necessary because 
     * gsl_eigen_symmv modifies it
    */
    gsl_matrix* m_working_mat;
    
    /*!
    * Used for storing the eigenvalues. We do not allocate memory for the pointer double*
    * because it will be directed to the right adress in the particle tag. The other 
    * information is set by hand in the constructor
    */
    gsl_vector m_evals;
    
    /*!
     * Used for storing the eigenvectors. We do not allocate memory for the pointer double*
     * because it will be directed to the right adress in the particle tag. The other 
     * information is set by hand in the constructor.
     */
    gsl_matrix m_evecs;
    
    /*!
    * The workspace for the computation of the eigensystem
    */
    gsl_eigen_symmv_workspace *m_eigen_workspace;
    
    /*!
     * Initialise the property list
     */
    virtual void init();

    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new PCaEigensystem(*this);
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
    
//  /*!
//   * Constructor
//   * @param color The particle's color
//   * @param offset Tag offset of the local volume
//   * @param tensor_offset Tag offset of the given tensor
    //   */
//    PCaEigensystem
//        (size_t color, size_t offset, string symbolName, size_t tensor_offset);
    
    /*!
     * Constructor for node hierarchy
     */
    PCaEigensystem(/*Node*/Simulation* parent);
 
  /*!
     * Destructor
   */
    virtual ~PCaEigensystem();

#ifdef _OPENMP
    /*!
     * Gets the number of doubles for a datatype, set by the \a DataFormat
     */
     virtual int setNumOfDoubles();
#endif 

  /*!
   * Compute the eigensystem
   */
   virtual void computeCacheFor(Particle* p) {/*
      p->tag.doubleByOffset(m_offset) += m_wf->interpolate(NULL, p->r) / p->tag.doubleByOffset(m_density_offset);*/
      
//         MSG_DEBUG("PCaEigensystem::computeCacheFor", "shear tensor: " << p->tag.tensorByOffset(8));
      
//         MSG_DEBUG("PCaEigensystem::computeCacheFor", "original tensor: " << p->tag.tensorByOffset(m_tensor_offset));
      
      // load the tensor
      tensor_t& tensor = p->tag.tensorByOffset(m_tensor_offset);
//       size_t k = 0;
//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "copied reference of tensor: " << tensor);
      
//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "GSL-working-tensor: ");
      for(size_t i = 0; i < SPACE_DIMS; ++i)
      {
        for(size_t j = 0; j < SPACE_DIMS; ++j)
        {
          GSL_MATRIX_SET(m_working_mat, i, j, tensor(i,j));
        
//           cout << GSL_MATRIX_GET(m_working_mat, i, j) << " ";
        }
//         cout << endl;
      } 
      
      // assign the right memory adresses to m_evals and m_evecs 
//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "evals in tag BEFORE: " << p->tag.pointByOffset(m_offset));
      
      m_evals.data = (double*)(p->tag.pointByOffset(m_offset)).coords;
      
/*      MSG_DEBUG("PCaEigensystem::computeCacheFor", "assigned m_evals BEFORE: ");
      for(size_t i = 0; i < SPACE_DIMS; ++i)
        cout << GSL_VECTOR_GET(&m_evals, i) << " ";
      cout << endl;
      MSG_DEBUG("PCaEigensystem::computeCacheFor", "evecs in tag BEFORE: " << p->tag.tensorByOffset(m_evecs_offset));*/
            
/*      p->tag.pointByOffset(m_offset)[0] = 1;
      p->tag.pointByOffset(m_offset)[1] = 2;
      p->tag.pointByOffset(m_offset)[2] = 3;
      MSG_DEBUG("PCaEigensystem::computeCacheFor", "evals in tag (point set for fun): " << p->tag.pointByOffset(m_offset));
      MSG_DEBUG("PCaEigensystem::computeCacheFor", "assigned m_evals (point set for fun): ");
      for(size_t i = 0; i < SPACE_DIMS; ++i)
        cout << GSL_VECTOR_GET(&m_evals, i) << " ";
      cout << endl;
      
      double* a = m_evals.data;
      *a = 10; ++a; *a = 20; ++a; *a = 30;
      MSG_DEBUG("PCaEigensystem::computeCacheFor", "evals in tag (m_evals set for fun): " << p->tag.pointByOffset(m_offset));
      MSG_DEBUG("PCaEigensystem::computeCacheFor", "assigned m_evals (m_evals set for fun): ");
      for(size_t i = 0; i < SPACE_DIMS; ++i)
        cout << GSL_VECTOR_GET(&m_evals, i) << " ";
      cout << endl;*/
 

      m_evecs.data =  ((p->tag.tensorByOffset(m_evecs_offset)).tensor);
      
/*      MSG_DEBUG("PCaEigensystem::computeCacheFor", "assigned m_evecs BEFORE: ");
      for(size_t i = 0; i < SPACE_DIMS; ++i)
      {
        for(size_t j = 0; j < SPACE_DIMS; ++j)
        {
          cout << GSL_MATRIX_GET(&m_evecs, i, j) << " ";
        }
        cout << endl; 
      }*/
      
      // Compute the eigensystem
      gsl_eigen_symmv(m_working_mat, &m_evals, &m_evecs, m_eigen_workspace);

//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "evals in tag AFTER: " << p->tag.pointByOffset(m_offset));
//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "assigned m_evals AFTER: ");
//       for(size_t i = 0; i < SPACE_DIMS; ++i)
//         cout << GSL_VECTOR_GET(&m_evals, i) << " ";
//       cout << endl;
//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "evecs in tag AFTER: " << p->tag.tensorByOffset(m_evecs_offset));
//       MSG_DEBUG("PCaEigensystem::computeCacheFor", "assigned m_evecs AFTER: ");
//       for(size_t i = 0; i < SPACE_DIMS; ++i)
//       {
//         for(size_t j = 0; j < SPACE_DIMS; ++j)
//         {
//           cout << GSL_MATRIX_GET(&m_evecs, i, j) << " ";
//         }
//         cout << endl; 
//       }

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
//       MSG_DEBUG("PCaEigensystem::==", "called");
      if (typeid(c) == typeid(*this)) {
        PCaEigensystem *cc = (PCaEigensystem*) &c;

        return
            m_offset == cc->m_offset && m_tensor_offset == cc->m_tensor_offset &&  m_tensor_symbol == cc->m_tensor_symbol &&  m_evecs_symbol == cc->m_evecs_symbol && m_evecs_offset == cc->m_evecs_offset && m_stage == cc->m_stage && m_offset == cc->m_offset && m_colour == cc->m_colour;
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
      temp.push_back(m_evecs_symbol);
      return temp;
    }

    /*!
     * Setup this Calculator
     */
    virtual void setup();
        
    /*!
     * Diffenrently to the function in \a Symbol, this class really has 
     * to determine its stage during run-time
     */
    virtual bool findStage();
        
    /*!
     * Diffenrently to the function in \a Symbol, this class really has 
     * to determine its stage during run-time
     */
    virtual bool findStage_0();
};

#endif
