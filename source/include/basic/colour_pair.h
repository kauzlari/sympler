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



#ifndef __COLOUR_PAIR_H
#define __COLOUR_PAIR_H 

#include "misc.h"
#include "gen_f.h"
#include "pair_list.h"
#include "data_format.h"
#include "manager_cell.h"

#ifdef _OPENMP
  #include "omp.h"
#endif


#define FOR_EACH_FREE_PAIR(cp, code)                                         \
{                                                                       \
  for (size_t t = 0; t < global::n_threads; ++t) {                  \
    SL_FOR_EACH(Pairdist, cp->freePairs()[t],\
		Pairdist* pair = __iSLFE;	   \
		code);			   \
  }						\
}		while(0)


#define FOR_EACH_PAIR(cp, code)                                         \
{                                                                       \
  for (int t = 0; t < global::n_threads; ++t) {                  \
if(cp->freePairsRandom(t).size())\
{\
  SL_FOR_EACH(PrimitiveSLEntry<size_t>, cp->freePairsRandom(t),\
      Pairdist* pair = &(cp->freePairs()[t][__iSLFE->m_val]);\
      code);\
}\
else \
{\
SL_FOR_EACH(Pairdist, cp->freePairs()[t],\
		Pairdist* pair = __iSLFE;	   \
		code);			   \
}\
if(cp->frozenPairsRandom(t).size()) \
{\
  SL_FOR_EACH(PrimitiveSLEntry<size_t>, cp->frozenPairsRandom(t),\
      Pairdist* pair = &(cp->frozenPairs()[t][__iSLFE->m_val]);\
      code);\
}\
else \
{\
  SL_FOR_EACH(Pairdist, cp->frozenPairs()[t],\
		Pairdist* pair = __iSLFE;	   \
		code);\
      }\
}						\
}		while(0)


#define FOR_EACH_PAIR_IN_THREAD(cp, t, code)                             \
{\
    SL_FOR_EACH(Pairdist, cp->freePairs(t),\
		Pairdist* pair = __iSLFE;	   \
		code);			   \
    SL_FOR_EACH(Pairdist, cp->frozenPairs(t),\
		Pairdist* pair = __iSLFE;	   \
		code);\
}		while(0)


#define FOR_EACH_PAIR_HALF_IN_GROUP(/*phase*/cp, groups, code)  \
{                                                               \
  for (int t = 0; t < global::n_threads; ++t) {              \
    if (groups.empty()) {                                       \
      SL_FOR_EACH(Pairdist, cp->freePairs()[t],			\
		  Pairdist* pair = __iSLFE;				\
		  double pair_factor = 1.0;			\
		  code);					\
    } else {                                                    \
      SL_FOR_EACH(Pairdist, cp->freePairs()[t],			\
		  Pairdist* pair = __iSLFE;				\
		  int gr1; int gr2;				\
		  						\
		  gr1 = pair->firstPart()->g;			\
		  gr2 = pair->secondPart()->g;			\
		  						\
		  bool a = groups.find(gr1) != groups.end();	\
		  bool b = groups.find(gr2) != groups.end();	\
		  if (a || b) {					\
		    double pair_factor = a && b ? 1.0 : 0.5;	\
		    code					\
		      }						\
		  );						\
    }								\
  }								\
} while(0)


/* ---- ColourPair ---- */

class Phase;

// class ManagerCell;


/*!
 * A \a ColourPair manages all properties of pairs of colors. This means
 * it stores, e.g., the list of all pairs of particles which have the color
 * combination this \a ColourPair refers to. A list of all \a ColourPair s
 * is managed by the \a ManagerCell. Note that color refers to a unique id of a
 * type of particle and species to the name of the type of the particle.
 */
class ColourPair
{
 protected:
  /*!
   * Pointer to the cell subdivision/pair list creation manager.
   */
  ManagerCell* m_manager;

  /*!
   * Does this \a ColourPair need to create pair lists? This might be
   * unnecessary, e.g., if this is a pair of two frozen colors.
   */
  bool m_needPairs;

  /*!
   * The \a ColourPair knows its place in the global ColourPair array.
   */
  size_t m_pos;

  /*!
   * The first color of this pair. Note that always \a m_1stColour < \a m_2ndColour.
   */
  size_t m_1stColour;

  /*!
   * The second color of this pair. Note that always \a m_1stColour < \a m_2ndColour.
   */
  size_t m_2ndColour;

  /*!
   * The cut-off radius for pair creation. This is set to the largest cut-off radius
   * that has been requested by any force, thermostat, etc.
   */
  double m_cutoff;

  /*!
   * Format description for the \a tag field
   */
  DataFormat m_tag_format;

  /*!
   * List of all pairs between free particles. The vector is needed such that
   * each thread has its own \a PairList so no concurrent writes to a \a PairList
   * are allowed to occur to avoid race conditions.
   */
  vector<PairList> m_freePairs;

  /*!
   * Indices of all pairs between free particles. Those maybe be ordered randomly 
   * for performing random loops over all pairs. The vector is needed such that
   * each thread has its own \a PairList so no concurrent writes to a \a PairList
   * are allowed to occur to avoid race conditions.
   */
  vector<SmartList<PrimitiveSLEntry<size_t> > > m_freePairsRandom;

  /*!
   * List of all pairs between free and frozen particles. The vector is needed such that
   * each thread has its own \a PairList so no concurrent writes to a \a PairList
   * are allowed to occur to avoid race conditions.
   */
  vector<PairList> m_frozenPairs;

  /*!
   * Indices of all pairs between free and frozen particles. Those maybe be ordered randomly 
   * for performing random loops over all pairs. The vector is needed such that
   * each thread has its own \a PairList so no concurrent writes to a \a PairList
   * are allowed to occur to avoid race conditions.
   */
  vector<SmartList<PrimitiveSLEntry<size_t> > > m_frozenPairsRandom;  

  /*!
   * Lists of connected pairs. The vector is for the different lists. A second one would be required for parallelisation
   */
  vector<PairList*> m_connectedLists;

  /*!
   * String identifiers of the connected lists in \a m_connectedLists 
   */
  vector<string> m_connectedListNames;
  
  /*!
   * A list of all forces that act between these two colors.
   */
  vector<GenF*> m_pair_forces;

  /*!
   * A list of all \a ValCalculator s for this color combination. A calculator
   * pre-calculates quantities that are used more than once, e.g., the
   * free energy or local density per particle which is needed more
   * than once for each particle. The different stages (second array)
   * are required so that
   * quantities that rely on other calculator quantities, which are only
   * available after having looped over all pairs once, can be calculated.
   * This does for example apply to the local shear rate which relies
   * on the local density.
   */
  vector<vector<ValCalculator*> > m_valCalculators;

  /*!
   * See documentation for \a m_valCalculators. These Calculators are for 
   * bonded pairs. An additional vector-layer is used for distinguishing 
   * between the different bonded lists available for one \a ColourPair
   */
  vector<vector<vector<ValCalculator*> > > m_bondedValCalculators;

  /*!
   * As \a m_bondedValCalculators but called at another instant during one timestep by the \a Controller
  */
  vector<vector<vector<ValCalculator*> > > m_bondedValCalculators_0;


#ifdef _OPENMP
  /*!
   * A list of all ValCalculatorPart's. for this colour combination.
   */
  vector<vector<ValCalculator*> > m_valCalculatorParts;
#endif

//  vector<ValCalculator*> m_valCalculators[/*N_CALCULATOR_STAGES*/VC_MAX_STAGE+1];

  /*!
   * This is a helper for building \a m_valCalculators . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_valCalculators .
   */
  vector<ValCalculator*> m_valCalculators_flat;

  /*!
   * This is a helper for building \a m_bondedValCalculators . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_bondedValCalculators .
   */
  /*vector<*/vector<ValCalculator*> /*>*/ m_bondedValCalculators_flat;

  /*!
   * The biggest stage occuring in \a m_valCalculators . This is determined during runtime in ColourPair::sortStages()
   */
  size_t m_maxStage;

  /*!
   * The biggest stage occuring in \a m_bondedValCalculators . This is determined during runtime in ColourPair::sortStages()
   */
  vector<size_t> m_maxBondedStage;

  /*!
   * As \a m_valCalculators but called at another instant during one timestep by the \a Controller
  */
  vector<vector<ValCalculator*> > m_valCalculators_0;

#ifdef _OPENMP
  /*!
   * As \a m_valCalculatorParts but called at another instant during one timestep by the \a Controller
   */
  vector<vector<ValCalculator*> > m_valCalculatorParts_0;
#endif

//  vector<ValCalculator*> m_valCalculators[/*N_CALCULATOR_STAGES*/VC_MAX_STAGE+1];

  /*!
   * This is a helper for building \a m_valCalculators_0 . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_valCalculators_0 .
   */
  vector<ValCalculator*> m_valCalculators_flat_0;

  /*!
   * This is a helper for building \a m_bondedValCalculators_0 . First, the \a ValCalculator s
   * are registered here and afterwards sorted by stages into \a m_bondedValCalculators_0 .
   */
  /*vector<*/vector<ValCalculator*> /*>*/ m_bondedValCalculators_flat_0;

  /*!
   * The biggest stage occuring in \a m_valCalculators_0 for one of the existing \a ColourPair s. This is determined during runtime in ColourPair::sortStages_0()
   */
  size_t m_maxStage_0;

  /*!
   * The biggest stage occuring in \a m_bondedValCalculators_0 . This is determined during runtime in ColourPair::sortStages_0()
   */
  vector<size_t> m_maxBondedStage_0;

 public:
  
   /*!
   * Standard constructor. Cannot be used.
   */
  ColourPair() {
    throw gError("ColourPair::ColourPair", "Don't call me! Contact the programmer.");
  }
  
  /*!
   * Constructor
   * @param manager Pointer to the cell subdivision manager
   */
  ColourPair(ManagerCell* manager)
    :  m_manager(manager), m_needPairs(false), m_1stColour(0), m_2ndColour(0), 
    m_cutoff(0),
    m_freePairs(global::n_threads, PairList(this)),
    m_freePairsRandom(global::n_threads, SmartList<PrimitiveSLEntry<size_t> >()),
    m_frozenPairs(global::n_threads, PairList(this)),
    m_frozenPairsRandom(global::n_threads, SmartList<PrimitiveSLEntry<size_t> >()) {
    }

  /*!
   * Copy constructor
   * @param cp \a ColourPair to copy
   */
  ColourPair(const ColourPair& cp);

  /*!
   * The position of this \a ColourPair in the global array of ColourPair's.
   */
  size_t &posInArr() {
    return m_pos;
  }

  vector<GenF*> *pairForces() {
    return &m_pair_forces;
  }

  /*!
   * return the \a maxStage 
   */
  size_t maxStage() const
  { 
    return m_maxStage;
  }
      
  /*!
   * return the \a maxStage 
   */
  size_t maxStage_0() const
  { 
    return m_maxStage_0;
  }
      
  /*!
   * return the \a m_maxBondedStage[icl] 
   */
  size_t maxBondedStage(size_t icl) const
  { 
    return m_maxBondedStage[icl];
  }
      
  /*!
   * return the \a m_maxBondedStage_0[icl] 
   */
  size_t maxBondedStage_0(size_t icl) const
  { 
    return m_maxBondedStage_0[icl];
  }
      
			
  /*!
   * Register a calculator in \a m_valCalculators_flat that caches quantities during pair creation.
   * These calculator store information in the \a Particle.
   * @param theSlots Return the tag offset of the quantities in the
   * respective \a Particle
   * @param vC The calculator
   * @param oneProp Is the calculated property specific for this ColourPair 
   * only?
   */ 
  void registerCalc(pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp);

  /*!
   * Register a calculator in \a m_valCalculators_flat_0 that caches quantities during pair creation.
   * These calculator store information in the \a Particle.
   * @param theSlots Return the tag offset of the quantities in the
   * respective \a Particle
   * @param vC The calculator
   * @param oneProp Is the calculated property specific for this ColourPair 
   * only?
   */ 
  void registerCalc_0(pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp);

			
  /*!
   * Register a calculator in \a m_bondedValCalculators_flat that caches quantities during pair creation.
   * These calculator store information in the \a Particle.
   * @param icl Index of the connected list of this calculator
   * @param theSlots Return the tag offset of the quantities in the
   * respective \a Particle
   * @param vC The calculator
   * @param oneProp Is the calculated property specific for this ColourPair 
   * only?
   */ 
  void registerBondedCalc(/*size_t icl,*/ pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp);

  /*!
   * Register a calculator in \a m_bondedValCalculators_flat_0 that caches quantities during pair creation.
   * These calculator store information in the \a Particle.
   * @param icl Index of the connected list of this calculator
   * @param theSlots Return the tag offset of the quantities in the
   * respective \a Particle
   * @param vC The calculator
   * @param oneProp Is the calculated property specific for this ColourPair 
   * only?
   */ 
  void registerBondedCalc_0(/*size_t icl,*/ pair<size_t, size_t> &theSlots, ValCalculator* vC, bool oneProp);

  /*!
   * Register a calculator in \a m_valCalculators_flat that caches quantities during pair creation.
   * These calculator store information in the \a Pardist. This only one offset
   * needs to be returned.
   * @param slot Return the tag offset of the quantities in the respective \a Pairdist
   * @param vC The calculator
   * @param oneProp Is the calculated property specific for this ColourPair 
   * only?
   */ 
  void registerCalc(size_t& slot, ValCalculator* vC, bool oneProp);
  
  /*!
   * Register a calculator in \a m_valCalculators_flat_0 that caches quantities during pair creation.
   * These calculator store information in the \a Pardist. This only one offset
   * needs to be returned.
   * @param slot Return the tag offset of the quantities in the respective \a Pairdist
   * @param vC The calculator
   * @param oneProp Is the calculated property specific for this ColourPair 
   * only?
   */ 
  void registerCalc_0(size_t& slot, ValCalculator* vC, bool oneProp);

  
  /*!
   * Register a calculator in \a m_valCalculators_flat that caches quantities during pair creation.
   * These calculator store information in the \a Pardist. This only one offset
   * needs to be returned.
   * This function assumes that the Valcalculator does all the registering work 
   * concerning the atribute itself and that the \a ValCalculator shoud definitely 
   * be registered in the ColourPair
   * @param vC The calculator
   */ 
  void registerCalc(ValCalculator* vC);
  	
  /*!
   * Register a calculator in \a m_valCalculators_flat_0 that caches quantities during pair creation.
   * These calculator store information in the \a Pardist. This only one offset
   * needs to be returned.
   * This function assumes that the Valcalculator does all the registering work 
   * concerning the atribute itself and that the \a ValCalculator shoud definitely 
   * be registered in the ColourPair
   * @param vC The calculator
   */ 
  void registerCalc_0(ValCalculator* vC);
  
  /*!
   * Register a calculator in \a m_bondedValCalculators_flat that caches quantities.
   * These calculator store information in the \a Pardist.
   * This function assumes that the Valcalculator does all the registering work 
   * concerning the atribute itself and that the \a ValCalculator shoud definitely 
   * be registered in the ColourPair
   * @param vC The calculator
   */ 
  void registerBondedCalc(ValCalculator* vC);
  	
  /*!
   * Register a calculator in \a m_bondedValCalculators_flat_0 that caches quantities.
   * These calculators store information in the \a Pardist.
   * This function assumes that the Valcalculator does all the registering work 
   * concerning the atribute itself and that the \a ValCalculator shoud definitely 
   * be registered in the ColourPair
   * @param vC The calculator
   */ 
  void registerBondedCalc_0(ValCalculator* vC);
  	
  /*!
   * Let the \a ValCalculator s in \a m_vlCalculators_flat determine their stages
   */
  bool findStages();
  
  /*!
   * Let the \a ValCalculator s in \a m_vlCalculators_flat_0 determine their stages
   */
  bool findStages_0();
  
  /*!
   * Sort the \a ValCalculator s from \a m_vlCalculators_flat into 
   * \a m_valCalculators according to their stages
   */
  void sortStages();
  
  /*!
   * Sort the \a ValCalculator s from \a m_vlCalculators_flat_0 into 
   * \a m_valCalculators_0 according to their stages
   */
  void sortStages_0();
  
  /*!
   * Return the number of calculators registered for stage \a stage.
   * @param stage Calculator stage
   */
  size_t nCalculatorsForStage(size_t stage) {
    return m_valCalculators[stage].size();
  }

  /*!
   * Register force \a force with this \a ColourPair. Usually invoked by
   * the force itself.
   * @param force Force to register
   */
  void registerForce(GenF* force) {
    m_pair_forces.push_back(force);
  }

  /*!
   * Return a pointer to the cell subdivision manager this \a ColourPair belongs to.
   */
  ManagerCell *manager() {
    return m_manager;
  }

  /*!
   * Does this \a ColourPair need pair calculation?
   */
  bool needPairs() const {
    return m_needPairs;
  }

  /*!
   * Set whether this \a ColourPair needs pair calculation
   */
  void setNeedPairs(bool np) {
    if (np)
      m_needPairs = true;
  }

  /*!
   * Return the first color of this pair.
   */
  size_t firstColour() const {
    return m_1stColour;
  }

  /*!
   * Return the second color of this pair.
   */
  size_t secondColour() const {
    return m_2ndColour;
  }

  /*!
   * Return the first species, i.e., the name and not the id, of this pair.
   */
  const string &firstSpecies() const;

  /*!
   * Return the second species, i.e., the name and not the id, of this pair.
   */
  const string &secondSpecies() const;

  /*!
   * Return the cut-off radius.
   */
  double cutoff() const {
    return m_cutoff;
  }

  /*!
   * Set the cut-off radius. Note: The cut-off radius cannot be decreased.
   */
  void setCutoff(double c);

  /*!
   * Returns a reference to the \a DataFormat of the \a tag.
   */
  DataFormat &tagFormat() {
    return m_tag_format;
  }
  
  /*!
   * Returns a reference to the list of free pairs for thread \a t
   * @param t Thread number
   */
  vector<PairList> &freePairs() {
    return m_freePairs;
  }
  
  /*!
   * Returns a reference to the list of frozen pairs for thread \a t
   * @param t Thread number
   */
  vector<PairList> &frozenPairs() {
    return m_frozenPairs;
  }

  /*!
   * Returns a reference to the ramdomised list of free pair indices for thread \a t
   * @param t Thread number
   */
  SmartList<PrimitiveSLEntry<size_t> > &freePairsRandom(size_t t) {
    return m_freePairsRandom[t];
  }
  
  /*!
   * Returns a reference to the randomised list of frozen pair indices for thread \a t
   * @param t Thread number
   */
  SmartList<PrimitiveSLEntry<size_t> > &frozenPairsRandom(size_t t) {
    return m_frozenPairsRandom[t];
  }

  /*!
   * Returns a reference to the list of frozen or free pairs for id \a t.
   * This means, the thread number is t mod number_of_thread, and if t/number_of_threads
   * is zero the list of free pairs is returned, otherwise the list of frozen pairs.
   * This is used to easily loop over all pairs, whether free or frozen. Fixme!!!
   * Perhaps get rid of frozen pairs alltogether and just disable 
   * the \a IntegratorPosition?
   */
/*   PairList &pairs(size_t t) { */
/*     if (t < global::n_threads) */
/*       return m_freePairs[t]; */
/*     else */
/*       return m_frozenPairs[t-global::n_threads]; */
/*   } */
  
  /*!
   * Return the list of calculators with index \a listIndex 
   * from \a m_bondedValCalculators for stage \a stage. 
   */
  vector<ValCalculator*>& bondedValCalculators(size_t stage, size_t listIndex) {
    return (m_bondedValCalculators[listIndex][stage]);
  }

  /*!
   * Return the list of calculators with index \a listIndex 
   * from \a m_bondedValCalculators_0 for stage \a stage. 
   */
  vector<ValCalculator*>& bondedValCalculators_0(size_t stage, size_t listIndex) {
    return (m_bondedValCalculators_0[listIndex][stage]);
  }

  /*!
   * Return the list of calculators from \a m_valCalculators for stage \a stage. 
   */
  vector<ValCalculator*> &valCalculators(size_t stage = 0) {
    return m_valCalculators[stage];
  }

  /*!
   * Return the list of calculators from \a m_valCalculators_0 for stage \a stage. 
   */
  vector<ValCalculator*> &valCalculators_0(size_t stage = 0) {
    return m_valCalculators_0[stage];
  }

#ifdef _OPENMP
  /*!
   * Return the list of calculators from \a m_valCalculatorParts for stage \a stage. 
   */
  vector<ValCalculator*> &valCalculatorParts(size_t stage = 0) {
    return m_valCalculatorParts[stage];
  }

  /*!
   * Return the list of calculators from \a m_valCalculatorParts_0 for stage \a stage. 
   */
  vector<ValCalculator*> &valCalculatorParts_0(size_t stage = 0) {
    return m_valCalculatorParts_0[stage];
  }
#endif

  /*!
   * Return the list of all calculators from \a m_valCalculators_flat
   */
  vector<ValCalculator*> &valCalculatorsFlat() {
    return m_valCalculators_flat;
  }

  /*!
   * Return the list of all calculators from \a m_valCalculators_flat_0
   */
  vector<ValCalculator*> &valCalculatorsFlat_0() {
    return m_valCalculators_flat_0;
  }

  /*!
   * Return the list of all calculators from \a m_bondedValCalculators_flat
   */
  vector<ValCalculator*>& bondedValCalculatorsFlat() {
    return m_bondedValCalculators_flat;
  }

  /*!
   * Return the list of all calculators from \a m_bondedValCalculators_flat_0
   */
  vector<ValCalculator*>& bondedValCalculatorsFlat_0() {
    return m_bondedValCalculators_flat_0;
  }

  /*!
   * Return a string identifier for this colour pair.
   */
  string toString() const {
    return "_" + firstSpecies() + "_" + secondSpecies() + "_";
  }

  /*!
   * Return a pointer to \a m_connectedLists
   */
  vector<PairList*>* connectedLists() {
    return &m_connectedLists;
  }

  /*!
   * Return the connected list from \a m_connectedLists at slot \a slot.
   */
  PairList* connectedList(size_t slot) {
    return m_connectedLists[slot];
  }
  
  /*!
   * Return the connected list from \a m_connectedLists corresponding to the identifier \a name.
   */
  PairList* connectedList(string name);
  
  /*!
   * Return the index of the connected list from \a m_connectedLists corresponding to the identifier \a name.
   */
  size_t connectedListIndex(string name);

  /*!
   * Return the index of the connected list from \a m_connectedLists corresponding to the identifier \a name. If the list does not exist yet it is created and written to a file.
   */
  size_t createConnectedListIndexAndWrite(string name);

  /*!
   * Return the index of the connected list from \a m_connectedLists corresponding to the identifier \a name. If the list does not exist yet it is created.
   */
  size_t createConnectedListIndex(string name);
  
  int searchConnectedListIndex(string name);
  
  /*!
   * Return the identifier of the connected list from \a m_connectedLists corresponding to the index \a name.
   */
  string connectedListName(size_t listIndex);
  
  /*!
   * Add the two particles to the connected list given by \a listIndex and write into a file
   */
  void addPairToConnectionAndWrite(Particle* p1, Particle* p2, size_t listIndex);

  /*!
   * Add the two particles to the connected list given by \a listIndex
   */
  void addPairToConnection(Particle* p1, Particle* p2, size_t listIndex);

  /*!
   * Does what its name indicates
   */
  void updateConnectedDistances();


  friend class ManagerCell;
//   friend class PairCreator;
};

#endif
