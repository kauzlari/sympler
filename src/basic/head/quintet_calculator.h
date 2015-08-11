/*
 * This file is part of the SYMPLER package.
 * https://github.com/kauzlari/sympler
 *
 * Copyright 2002-2015, 
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


#ifndef __QUINTET_CALCULATOR_H
#define __QUINTET_CALCULATOR_H

#include "symbol.h"

class Pairdist;
class ColourPair;

/*!
 * Parent class for a \a Symbol for cached properties computed during a loop over bonded quintets
 * NOTICE: This Calculator relies on convention 6 in the CONVENTIONS file
 */
class QuintetCalculator : public Symbol
{
protected:
  
  /*!
   * Should periodic boundary conditions be applied to the connections?
   */
  bool m_periodic;

  /*! helper: boxSize for determining periodic boundary conditions*/
  point_t m_boxSize;

  /*!
   * index of the bonded list, this calculator belongs to
   */
  size_t m_listIndex;
    
  /*!
   * name of the bonded list, this calculator belongs to
   */
  string m_listName;

  /*!
   * Memory offsets to the symbol in the \a Particle 's tag this calculator computes, one for each of the possibly different species 
   */
  size_t m_slots[5];

#ifdef _OPENMP
  /*!
   * The tag offset of the copy data stored in a \a Particle.
   * The vector size equals the number of threads used. 
   * The size_t* should point to arrays of three size_t
   */
    vector<size_t*> m_copy_slots;

    /*!
     * The offset to the copy-slots inside a vector
     */
    size_t m_vector_slots[5];
#endif

  /*!
   * Species of the quintet particles. This relies on CONVENTION 6.
   */
  string m_species[5];

  /*!
   * Colour of first \a Particle in the Quintet
   */
  size_t m_firstColour;

  /*!
   * Colour of first \a Particle in the Quintet
   */
  size_t m_secondColour;

  /*!
   * Colour of first \a Particle in the Quintet
   */
  size_t m_thirdColour;

  /*!
   * Colour of first \a Particle in the Quintet
   */
  size_t m_fourthColour;

  /*!
   * Colour of first \a Particle in the Quintet
   */
  size_t m_fifthColour;

    /*! This is now a member of \a Symbol
     * Is this calculator allowed to overwrite already existing symbols 
     * with name \a m_symbolName ?
     */
/*     bool m_overwrite; */

#ifdef _OPENMP
  /*!
   * Is this a particle or a pair calculator.
   */
  bool m_particleCalculator;
#endif
  
  /*!
  * Initialise the property list
  */
  virtual void init();

  /*!
   * Helper function for polymorphic copying
   */
  virtual QuintetCalculator* copyMySelf() = 0;
/*   { */
/*     return new TripletCalculator(*this); */
/*   } */


public:
  /*!
 * Constructor
   */
  QuintetCalculator(string symbol);
  
  /*!
   * Constructor
   */
  QuintetCalculator(/*Node*/Simulation* parent);
   
  /*!
   * Destructor
   */
  virtual ~QuintetCalculator() {
  }

  /*!
   * Compute cached properties for \a triplet_t \a tr
   * @param tr The \a triplet_t for which to compute the cached properties
   */
#ifndef _OPENMP
  virtual void compute(quintet_t* tr) = 0;
#else
  virtual void compute(quintet_t* tr/*, size_t thread_no*/ /*FIXME: parallelise!*/) = 0;
#endif

#ifdef _OPENMP
  virtual bool particleCalculator() {
    return m_particleCalculator; 
  }
#endif


  /*!
   * Return a string identifier for this calculator
   */
  virtual string myName() {
    return m_symbolName;
  }

#ifdef _OPENMP
  /*!
   * Return the offset for the data stored by this calculator in the \a Particle s
   */
    virtual vector<size_t*> &copySlots() {
      return m_copy_slots;
    }

  /*!
   * Return the offset for the data inside the copy vector in the \a Particle s tag.
   */
    virtual size_t* vectorSlots() {
      return m_vector_slots;
    }

  /*!
   * Merging the copies among different threads (processors) together
   */
    virtual void mergeCopies(size_t thread_no) = 0;
#endif

    /*!
     * setup this calculator
     */
  virtual void setup();

    /*!
     * run those setups requiring that \a Particle s have been created
     */
  virtual void setupAfterParticleCreation();

  /*!
   * Returns the name of the bonded list this calculator computes for.
   */
  virtual string listName() {
    return m_listName;
  }

    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
  virtual bool findStage();
    
    /*!
     * Determines \a m_stage of the current \a Symbol.
     * By default, we assume that the stage is fixed and known during compile-time, 
     * so this function does nothing except returning the message (true) that the 
     * stage was already found. Symbols, which determine the stage during run-time 
     * have to redefine this function.
     */
    virtual bool findStage_0();

    /*!
     * Checks for a consistent setup. Currently (2009-08-10) 
     * just needed for the quintetCalculator 
     */
    virtual void checkConsistency();

};


#endif
