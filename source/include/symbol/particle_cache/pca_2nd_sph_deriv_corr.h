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


#ifndef __PCA_2ND_DERIV_CORR_H
#define __PCA_2ND_DERIV_CORR_H


#include "gsl_helper.h"
#include "weighting_function.h"

#include "particle_cache.h"


/*!
 * Computes the symmetric 3x3 correction matrix B for the second SPH-derivative
 * according to R. Fatehi, M.T. Manzari, Computers and Mathematics with
 * Applications 61 (2011) 482–498. B is obtained from the equation B:A == -I,
 * where A is a 3x3x3x3 matrix with some symmetries, which allows A to be
 * represented as a non-symmetric 6x6-matrix, while B is stored in a 6-vector.
 * I is the 3x3 identity matrix.
 */
class PCa2ndSPHDerivCorr: public ParticleCache
{

  protected:

	  /*!
	   * The \a WeightingFunction to be used. This class will use the
	   * WeightingFunction::weight() method
	   */
	  WeightingFunction *m_wf;

	  /*!
	   * Name of the \a WeightingFunction \a m_wf
	   */
	  string m_wfName;

		/*!
		 * Tag offset of the externally computed particle volume. This is a vector
		 * because we need one per species we are interacting with. Storage logic:
		 * First entry: own species. Following entries: other species in the order
		 * we access them when looping over \a ColourPair s.
		 */
		vector<size_t> m_volumeOffset;

    /*!
     * Symbol name of the externally computed particle volume
     */
    string m_volumeName;

		/*!
		 * Tag offset of the internally computed symmetric 3x3x3x3 matrix A of the
		 * linear system B:A == I defining the correction matrix B. Due to its
		 * symmetry, A is treated as a non-symmetric 6x6 matrix and stored
		 * row-major in a vector<double>.
		 * NOTE: no string required for storing internal hard-coded symbol name.
		 * FIXME: this means we store 36 doubles simultaneously for each
		 * \a Particle!! This is a horrible waste of memory since we could also
		 * solve the equation system immediately and reuse the same memory for the
		 * next particle! Reason for the waste: The current (2018-07-31)
		 * \a PairList does not permit for a true neighbour list for each
		 * \a Particle, i.e., to loop over all the \a Particle's neighbours, and
		 * only those, quickly!
		 */
		size_t m_systemMatOffset;

		/*!
		 * Tag offset of the externally computed correction matrix for the
		 * SPH-discretisation of the 1st derivative
		 */
		size_t m_1stDerivCorrOffset;

		/*!
		 * Symbol name of the externally computed correction matrix for the
		 * SPH-discretisation of the 1st derivative
		 */
		string m_1stDerivCorrName;

		/*!
		 * Tag offset for a 6x6 helper matrix Cabcd. See code for its purpose.
		 * NOTE: no string required for storing internal hard-coded symbol name.
		 */
		size_t m_CabcdHelperOffset;

		/*!
		 * Tag offset for a 6x3 helper matrix Habf. See code for its purpose.
		 * NOTE: no string required for storing internal hard-coded symbol name.
		 */
		size_t m_HabfHelperOffset;

		/*!
		 * Tag offset for a 6x6 helper matrix Eabe. See code for its purpose.
		 * NOTE: no string required for storing internal hard-coded symbol name.
		 */
		size_t m_EabeHelperOffset;

		/*!
		 * Tag offset for a 6x6 helper matrix T_Gcab_abc. See code for its purpose.
		 * NOTE: no string required for storing internal hard-coded symbol name.
		 */
		size_t m_T_Gcab_abcHelperOffset;

    /*!
     * The workspace for the computation of the LU decomposition
     */
    gsl_permutation *m_permutation;
    
    /*!
     * GSL handle of the equation system matrix B. No dynamic allocation
     * needed, since we can simply assign memory addresses to the previously
     * computed matrix pointed to by \a m_systemMatOffset (using the current
     * bad technique, FIXME!), and since it is no problem that the matrix gets
     * modified during LU decomposition. The other required information is set
     * by hand in the constructor.
     */
    gsl_matrix m_workingMat;
    
    /*!
     * Storage for the solution B of B:A == -I in 6-vector form. Will be copied
     * afterwards into the symmetric 3x3 tensor stored in \a Particle::tag at
     * memory offset \a m_offset.
     * This seems to be the price to pay if we want to reduce the system from
     * 3x3x3x3 = 9x9 to 6x6 due to symmetry.
     * FIXME: Is it worth it? Any better solution? Well, we could introduce a
     * DataFormat::SYMTENSOR holding a double[6] array, including a way to use
     * them in runtime compiled expressions in the same way as standard TENSORs
     * Storage convention: for B-vector in terms of 3x3 B-matrix Bab:
     * B = (B11, B22, B33, 2*B12, 2*B13, 2*B23)
     * The 3 factors of 2 here (+3 divisions by 2 later) avoid 18 of these
     * factors in the A-matrix.
     */
    gsl_vector* m_B;
    
    /*!
     * Helper used to compute once the A-matrices for all particles by pair
     * summation at the first call of computeCacheFor for the first particle.
     * This boolean is set to true in precompute() and to false in
     * computeCacheFor().
     */
    bool m_pairLoopToDo;

    /*!
     * The right-hand-side as 6-vector, with 3 zeros and 3 (minus) ones
     */
    static gsl_vector* s_rhs;

    /*!
     * Number of entries of the reduced B matrix (6 in 3D)
     * = SPACE_DIMS*(SPACE_DIMS+1)/2
     */
    static const size_t s_matEntries;

    /*!
     * Initialise the property list
     */
    virtual void init();

    /*!
     * Helper function for polymorphic copying
     */
    virtual ParticleCache* copyMySelf()
    {
      return new PCa2ndSPHDerivCorr(*this);
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
      usedSymbols.push_back(m_volumeName);
      usedSymbols.push_back(m_1stDerivCorrName);
    }

    /*!
     * Computes the matrix A of the linear system B . A = -I for each particle
     */
		virtual void buildSystemMatrices();

    /*!
     * Initialises all matrices including the helper matrices to zero
     */
		virtual void initSystemMatrices();


  public:

    /*!
     * Constructor for node hierarchy
     */
    PCa2ndSPHDerivCorr(/*Node*/Simulation* parent);

    /*!
     * Destructor
     */
    virtual ~PCa2ndSPHDerivCorr();

    /*!
     * Compute the correction matrix
     */
    virtual void computeCacheFor(Particle* p);

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
        PCa2ndSPHDerivCorr *cc = (PCa2ndSPHDerivCorr*) &c;

        return
        		// FIXME: this is ridiculous! Do we still need this non-sense?
            m_offset == cc->m_offset && m_wf == cc->m_wf &&
						m_wfName == cc->m_wfName && m_volumeOffset == cc->m_volumeOffset &&
						m_volumeName == cc->m_volumeName &&
						m_systemMatOffset == cc->m_systemMatOffset &&
						m_1stDerivCorrOffset == cc->m_1stDerivCorrOffset &&
						m_CabcdHelperOffset == cc->m_CabcdHelperOffset &&
						m_HabfHelperOffset == cc->m_HabfHelperOffset &&
						m_EabeHelperOffset == cc->m_EabeHelperOffset &&
						m_T_Gcab_abcHelperOffset == cc->m_T_Gcab_abcHelperOffset &&
						m_permutation == cc->m_permutation &&
						&m_workingMat == &(cc->m_workingMat) &&
						m_B == cc->m_B &&
						m_pairLoopToDo == cc->m_pairLoopToDo;
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

    /*!
     * Checks existence of input symbols required by this \a ParticleCache
     * in a hard-coded fashion (i.e., not through runtime compiled
     * expressions).
     * @param colour Particle colour to be checked
     */
    virtual void checkInputSymbolExistences(size_t colour);

    /*!
     * Runs required operations before the actual computations in
     * computeCacheFor()
     */
    virtual void precompute();

    /*!
     * Return static member \a s_matEntries
     */
    static size_t return_s_matEntries() {

    	return s_matEntries;
    }

    /*!
     * Return static member \a s_rhs
     */
    static gsl_vector* return_s_rhs() {

    	return s_rhs;
    }

    /*!
     * Nested class for initialising static members of class
     * \a PCa2ndSPHDerivCorr since some of the initialisation requires
     * non-trivial code not writable as a one-line.
     */
    class _Initialiser
		{
    	public:

    		/*!
    		 * The default constructor will initialize our static variable
    		 */
    		_Initialiser()
 				{
    			for (size_t i = 0; i < PCa2ndSPHDerivCorr::s_matEntries; ++i) {

    				double val;
    				if (i < SPACE_DIMS) val = -1.0;
    				else val = 0.0;

    				gsl_vector_set(PCa2ndSPHDerivCorr::s_rhs, i, val);

    			}
 				}
		};

	private:

    /*!
     * We'll use this static object to ensure the
     * \a _Initialiser::_Initialiser() constructor is called
     */
    static _Initialiser s_initialiser;

    /*!
     * Helper storage for matrix C of "my" particle for usage across functions
     */
//    vector<double>& myCabcd;
    /*!
     * Helper storage for matrix H of "my" particle for usage across functions
     */
//    vector<double>& myHabf;
    /*!
     * Helper storage for matrix E of "my" particle for usage across functions
     */
//    vector<double>& myEabe;
    /*!
     * Helper storage for matrix G of "my" particle for usage across functions
     */
//    vector<double>& myT_Gcab_abc;
    /*!
     * Helper storage for matrix C of "other" particle for usage across
     * functions
     */
//    vector<double>& otherCabcd;
    /*!
     * Helper storage for matrix H of "other" particle for usage across
     * functions
     */
//    vector<double>& otherHabf;
    /*!
     * Helper storage for matrix E of "other" particle for usage across
     * functions
     */
//    vector<double>& otherEabe;
    /*!
     * Helper storage for matrix G of "other" particle for usage across
     * functions
     */
//    vector<double>& otherT_Gcab_abc;
        
};

#endif


